use crate::simulation::accessgpuspins::{get_spin, use_offset};
use crate::simulation::interactions::GPUInteractions;
use crate::utils::constants::*;
use crate::utils::vec3::*;
use cubecl::prelude::*;
use crate::utils::gpurng::{gaussian_pair, rng_state};

#[cube]
fn multiply_tensor<F: Float>(
    interactions:&GPUInteractions,
    index:usize,
    spin:Vec3<F>
) -> Vec3<F> {
    let jxx = F::cast_from(interactions.jxx[index]);
    let jxy = F::cast_from(interactions.jxy[index]);
    let jxz = F::cast_from(interactions.jxz[index]);
    let jyx = F::cast_from(interactions.jyx[index]);
    let jyy = F::cast_from(interactions.jyy[index]);
    let jyz = F::cast_from(interactions.jyz[index]);
    let jzx = F::cast_from(interactions.jzx[index]);
    let jzy = F::cast_from(interactions.jzy[index]);
    let jzz = F::cast_from(interactions.jzz[index]);

    Vec3::<F> {
        x: jxx*spin.x + jxy*spin.y + jxz*spin.z,
        y: jyx*spin.x + jyy*spin.y + jyz*spin.z,
        z: jzx*spin.x + jzy*spin.y + jzz*spin.z,
    }
}

#[cube]
pub fn get_interaction_field<F: Float>(
    interactions:&GPUInteractions,
    spins:&Array<F>,
    layer:u32,
    uc_x:u32,
    uc_y:u32,
    site:u32,
    layers:u32,
    x_cells:u32,
    y_cells:u32,
    site_magnetic_moment:f32,
    num_interactions:u32,
) -> Vec3<F> {
    let interaction_index = (layer as usize % 6)*4*num_interactions as usize 
        + (site * num_interactions) as usize;
    let mut bx = F::zero();
    let mut by = F::zero();
    let mut bz = F::zero();
    for pair_atom in 0..num_interactions as usize {
        let pair_spin = get_spin(
            spins,
            use_offset(layer, interactions.layer_offset[interaction_index+pair_atom], layers),
            use_offset(uc_x, interactions.uc_x_offset[interaction_index+pair_atom], x_cells),
            use_offset(uc_y, interactions.uc_y_offset[interaction_index+pair_atom], y_cells),
            interactions.pair_site[interaction_index+pair_atom], x_cells, y_cells);
        let tensor_result = multiply_tensor(
            interactions,
            interaction_index+pair_atom,
            pair_spin);
        let to_field = vec3_scalar_multiply(
            F::one()/F::cast_from(site_magnetic_moment*BOHR_MAGNETON),
            &tensor_result);
        bx += to_field.x;
        by += to_field.y;
        bz += to_field.z;
    }
    Vec3::<F> {
        x: bx,
        y: by,
        z: bz
    }
}

#[cube]
pub fn get_anisotropy_field<F: Float>(
    spin:&Vec3<F>,
    anisotropy_strength:F,
    site_magnetic_moment:f32,
) -> Vec3<F> {
    let anisotropy_axis = Vec3::<F> {
        x:F::zero(),
        y:F::zero(),
        z:F::one()
    };
    let prefactor = F::cast_from(2.0f32)* anisotropy_strength /F::cast_from(site_magnetic_moment*BOHR_MAGNETON);
    let dot = vec3_dot(spin, &anisotropy_axis);
    vec3_scalar_multiply(prefactor*dot, &anisotropy_axis)
}

#[cube]
pub fn get_external_field<F: Float>(
    uc_x:u32, uc_y:u32,
    x_cells: u32, y_cells: u32,
    phase_sin: F,
    sigma:F,
    oscillation_vector: &Vec3<F>,
    b_const: &Vec3<F>,
) -> Vec3<F>{
    let dx = (uc_x as i32 - (x_cells as i32 / 2)) as f32*5.105f32;
    let dy = (uc_y as i32 - (y_cells as i32 / 2)) as f32*8.842f32;
    let radius_sq = F::cast_from(dx * dx + dy * dy);
    let oscillate = phase_sin * (-radius_sq/(F::cast_from(2.0f32)*sigma*sigma)).exp();
    let b_osc = vec3_scalar_multiply(F::cast_from(oscillate), oscillation_vector);
    vec3_add(b_const, &b_osc)
}

#[cube]
pub fn get_thermal_field<F: Float>(
    layer:u32,
    uc_x:u32,
    uc_y:u32,
    site:u32,
    global_rng_state:u32,
    temperature_factor:F,
) -> Vec3<F> {
    let mut state = rng_state(layer, uc_x, uc_y, site);
    state ^= global_rng_state;
    let (g1,g2) = gaussian_pair(&mut state);
    let (g3, _g4) = gaussian_pair(&mut state);

    Vec3::<F> {
        x: F::cast_from(g1)*temperature_factor,
        y: F::cast_from(g2)*temperature_factor,
        z: F::cast_from(g3)*temperature_factor,
    }
}

#[cube]
pub fn get_spin_orbit_field<F: Float>(
    spin:&Vec3<F>,
    sigma:&Vec3<F>,
    b_damping_like:f32,
    b_field_like:f32
) -> Vec3<F> {
    let b_dl_vec = vec3_scalar_multiply(F::cast_from(b_damping_like), &vec3_cross(spin, sigma));
    let b_fl_vec = vec3_scalar_multiply(F::cast_from(b_field_like), sigma);

    vec3_add(&b_fl_vec, &b_dl_vec)
}