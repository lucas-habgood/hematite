use crate::simulation::accessgpuspins::{get_array_index, pos_from_index};
use crate::simulation::field::{get_external_field, get_thermal_field};
use crate::simulation::interactions::{GPUInteractions, PackedInteractions, interaction_launch, upload_interactions};
use crate::simulation::llg::llg_rk4;
use crate::statistics::logsamples::log_samples;
use crate::statistics::magvectors::get_mag_neel;
use crate::statistics::output::handle_output;
use crate::statistics::savestate::save_state_to_file;
use crate::statistics::topological::get_topological_charge;
use crate::utils::config::config;
use crate::utils::constants::{BOHR_MAGNETON, BOLTZMANN_CONSTANT, GYROMAGNETIC_RATIO};
use crate::utils::cpuspinarray::site_to_index;
use crate::utils::vec3::{Vec3, vec3_add};
use crate::visualisation::createplot::plot;
use cubecl::prelude::*;
use cubecl::{CubeElement, Runtime};
use rand::RngExt;
use rand::rngs::StdRng;
use std::f32::consts::TAU;

#[cube(launch_unchecked)]
fn step_spins<F: Float>(
    interactions_data: &GPUInteractions,
    spins: &Array<F>, out_spins: &mut Array<F>,
    mask: &Array<u32>,
    layers: u32, x_cells: u32, y_cells: u32,
    dispatch_width: u32,
    num_spins: u32,
    num_interactions: u32,
    timestep:f32,
    anisotropy_strength:f32,
    alpha:f32,
    applied_field_x:f32, applied_field_y:f32, applied_field_z:f32,
    site_magnetic_moment:f32,
    phase_sin:f32,
    osc_sigma:f32,
    osc_field_x:f32, osc_field_y:f32, osc_field_z:f32,
    global_rng_state:u32,
    temperature_factor:f32,
    sot_b_dl:f32, sot_b_fl:f32,
    sot_sigma_x:f32, sot_sigma_y:f32, sot_sigma_z:f32,
)
{
    let i = ABSOLUTE_POS_Y * dispatch_width + ABSOLUTE_POS_X;

    if i<num_spins && mask[i as usize]==1{
        let (layer, uc_x, uc_y, site) =
            pos_from_index(i, x_cells, y_cells);
        let b_applied = Vec3::<F> {
            x: F::cast_from(applied_field_x),
            y: F::cast_from(applied_field_y),
            z: F::cast_from(applied_field_z),
        };

        let oscillation_vector = Vec3::<F> {
            x: F::cast_from(osc_field_x),
            y: F::cast_from(osc_field_y),
            z: F::cast_from(osc_field_z),
        };

        let sot_sigma = Vec3::<F> {
            x: F::cast_from(sot_sigma_x),
            y: F::cast_from(sot_sigma_y),
            z: F::cast_from(sot_sigma_z),
        };
        // if this external field is expensive, and changes slowly relative to iteration timestep
        // then could just compute in another kernel every so often
        let b_ext = get_external_field(
            uc_x, uc_y,
            x_cells, y_cells,
            F::cast_from(phase_sin),
            F::cast_from(osc_sigma),
            &oscillation_vector,
            &b_applied,
        );

        let b_therm = get_thermal_field(
            layer, uc_x, uc_y, site,
            global_rng_state,
            F::cast_from(temperature_factor),
        );

        let updated_spin = llg_rk4(
            interactions_data, spins,
            layers, x_cells, y_cells,
            timestep, anisotropy_strength, alpha,
            layer, uc_x, uc_y, site,
            num_interactions,
            site_magnetic_moment,
            &vec3_add(&b_therm, &b_ext),
            sot_b_dl, sot_b_fl, &sot_sigma
        );

        let spindex = get_array_index(layer, uc_x, uc_y, site, x_cells, y_cells);
        out_spins[spindex] = updated_spin.x;
        out_spins[spindex + 1] = updated_spin.y;
        out_spins[spindex + 2] = updated_spin.z;
    }
}

#[cube(launch_unchecked)]
fn sample_spins<F: Float>(
    spins: &Array<F>,
    indices: &Array<u32>,
    sample: &mut Array<F>,
    write_index: u32,
) {
    let i = ABSOLUTE_POS;
    let index = indices[i] as usize;
    let sample_index = write_index as usize + i;

    sample[3 * sample_index] = spins[3 * index];
    sample[3 * sample_index + 1] = spins[3 * index + 1];
    sample[3 * sample_index + 2] = spins[3 * index + 2];

}

fn sample_integer_positions(
    samples: Vec<[f32;3]>,
    uc_dims: [usize; 3],
) -> Vec<[usize;3]> {
    let mut positions = Vec::<[usize;3]>::new();
    for sample in samples {
        let sample_uc_x = ((sample[0] * uc_dims[0] as f32).round() as usize).clamp(0, uc_dims[0]-1);
        let sample_uc_y = ((sample[1] * uc_dims[1] as f32).round() as usize).clamp(0, uc_dims[1]-1);
        let sample_layer = ((6.0*sample[2]*uc_dims[2] as f32).round() as usize).clamp(0, 6*uc_dims[2]-1);
        positions.push([sample_uc_x, sample_uc_y, sample_layer]);
    }
    positions
}


pub fn launch_gpu<R: Runtime>(
    uc_dims:[usize;3],
    interactions: PackedInteractions,
    spin_array: Vec<f32>,
    mask: Vec<u32>,
    rng: &mut StdRng,
) {


    let device = <R as Runtime>::Device::default();
    let client = R::client(&device);

    let x_cells = uc_dims[0] as u32;
    let y_cells = uc_dims[1] as u32;
    let layers = uc_dims[2] as u32 * 6u32;


    const MAX_DISPATCH_X: u32 = 65535;
    const WORKGROUP_SIZE: u32 = 256;

    let cube_dim = CubeDim::new_1d(WORKGROUP_SIZE);

    let num_spins = spin_array.len() as u32 / 3;
    let num_workgroups = num_spins.div_ceil(WORKGROUP_SIZE);

    let cube_count_x = num_workgroups.min(MAX_DISPATCH_X);
    let cube_count_y = num_workgroups.div_ceil(cube_count_x);

    let cube_count = CubeCount::Static(
        cube_count_x,
        cube_count_y,
        1,
    );

    let mut input_handle = client.create_from_slice(f32::as_bytes(&spin_array));
    let mut output_handle = client.empty(spin_array.len() * size_of::<f32>());
    let mask_handle = client.create_from_slice(u32::as_bytes(&mask));

    let sample_uc_coords = sample_integer_positions(
        config().samples.clone(), uc_dims
    );
    let mut sample_indices:Vec<u32> = vec![];
    for sample_uc_coord in sample_uc_coords {
        sample_indices.push(site_to_index(
            config().plot_layer, sample_uc_coord[0], sample_uc_coord[1], 0,
            x_cells as usize, y_cells as usize) as u32
        );
    }
    let num_samples = sample_indices.len();
    let indices_handle = client.create_from_slice(u32::as_bytes(&sample_indices));
    let sample_handle = client.empty( // batch the sample reads
          num_samples*3 * size_of::<f32>()*config().iters_per_plot as usize);

    unsafe {
        let interaction_storage = upload_interactions(
            &interactions, &client);
        let interactions_per_spin = interactions.layer_offset.len()/24;
        if config().precondition_iterations>0 {println!("Preconditioning");}
        for _pre_iter in 1..config().precondition_iterations+1 {
            let gpu_spin_array = ArrayArg::<R>::from_raw_parts(
                input_handle.clone(), spin_array.len());
            let output_array = ArrayArg::<R>::from_raw_parts(
                output_handle.clone(), spin_array.len());
            let gpu_interactions_data =
                interaction_launch(&interaction_storage, interactions_per_spin);
            let mask_array = ArrayArg::<R>::from_raw_parts(
                mask_handle.clone(), mask.len());

            step_spins::launch_unchecked::<f32, R>(
                &client, cube_count.clone(), cube_dim,
                gpu_interactions_data,
                gpu_spin_array,
                output_array,
                mask_array,
                layers, x_cells, y_cells,
                cube_count_x*WORKGROUP_SIZE, num_spins,
                interactions_per_spin as u32,
                config().precondition_timestep,
                config().anisotropy_strength,
                1.0, // high damping for precondition
                0.0, 0.0, 0.0, // dont apply any field while preconditioning
                config().site_magnetic_moment,
                0.0, // no oscillation
                1.0, // avoid divide by zero, no effect anyway
                0.0, 0.0, 0.0, // ^
                0, // thermal rng irrelevant
                0.0, // zero temperature
                0.0, 0.0, // no SOT
                0.0, 0.0, 0.0, // SOT direction irrelevant
            );
            std::mem::swap(&mut input_handle, &mut output_handle);
        }



        println!("Launching Kernel");
        let start = std::time::Instant::now();
        for iter in 1..config().iterations+1 {
            let temperature = config().initial_temperature +
                (config().final_temperature-config().initial_temperature) *
                    (iter as f32 / config().iterations as f32);
            let temperature_factor = (
                2.0*config().alpha*BOLTZMANN_CONSTANT*temperature
                    /(GYROMAGNETIC_RATIO*BOHR_MAGNETON*config().site_magnetic_moment*config().timestep)
            ).sqrt();

            let gpu_spin_array = ArrayArg::<R>::from_raw_parts(
                input_handle.clone(), spin_array.len());
            let output_array = ArrayArg::<R>::from_raw_parts(
                output_handle.clone(), spin_array.len());
            let gpu_interactions_data =
                interaction_launch(&interaction_storage, interactions_per_spin);
            let mask_array = ArrayArg::<R>::from_raw_parts(
                mask_handle.clone(), mask.len());

            step_spins::launch_unchecked::<f32, R>(
                &client, cube_count.clone(), cube_dim,
                gpu_interactions_data,
                gpu_spin_array,
                output_array,
                mask_array,
                layers, x_cells, y_cells,
                cube_count_x*WORKGROUP_SIZE, num_spins,
                interactions_per_spin as u32,
                config().timestep,
                config().anisotropy_strength,
                config().alpha,
                config().applied_field[0], config().applied_field[1], config().applied_field[2],
                config().site_magnetic_moment,
                (TAU*config().osc_freq *config().timestep*iter as f32).sin(),
                config().osc_sigma,
                config().osc_field[0], config().osc_field[1], config().osc_field[2],
                rng.random::<u32>(),
                temperature_factor,
                config().sot_b_dl, config().sot_b_fl, // SOT damping & field like effective B field
                config().sot_current_direction[0],
                config().sot_current_direction[1],
                config().sot_current_direction[2],
            );
            if num_samples>0 {
                sample_spins::launch_unchecked::<f32, R>(
                    &client, CubeCount::Static(1,1,1), CubeDim::new_1d(num_samples as u32),
                    ArrayArg::<R>::from_raw_parts(
                        output_handle.clone(), spin_array.len()),
                    ArrayArg::<R>::from_raw_parts(
                        indices_handle.clone(),
                        num_samples,
                    ),
                    ArrayArg::<R>::from_raw_parts(
                        sample_handle.clone(),
                        num_samples*3
                    ),
                    ((iter-1)%config().iters_per_plot) * num_samples as u32,
                );
            }

            if iter%config().iters_per_plot==0 {
                let bytes = client.read_one(output_handle.clone()).expect("Failed to read output");
                let output = f32::from_bytes(&bytes);
                let filename = format!("plots/{}.png", iter);
                if !config().stats_only {plot(&output, uc_dims, &*filename);}
                let (mag, neel) = get_mag_neel(&output, uc_dims, &mask);
                let charge = get_topological_charge(
                    &output, uc_dims, &mask
                );
                handle_output(iter, &start, temperature, charge, mag, neel).expect("Failed to output");

                if num_samples>0 {
                    let samples_bytes = client.read_one(sample_handle.clone())
                        .expect("Failed to read from sample handle");
                    let samples = f32::from_bytes(&samples_bytes);
                    log_samples(&*config().sample_output_file.clone(), samples,
                                config().iters_per_plot as usize, num_samples
                    ).expect("Failed to log samples");
                }
            }
            std::mem::swap(&mut input_handle, &mut output_handle);
        }

        // end save state
        if config().save_at_end {
            let bytes = client.read_one(output_handle.clone()).expect("Failed to read output");
            let output = f32::from_bytes(&bytes);
            save_state_to_file(
                output,
                uc_dims,
                &*config().save_file
            ).expect("Failed to save state");
        }
    }
}