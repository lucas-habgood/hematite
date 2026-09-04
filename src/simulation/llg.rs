use crate::utils::constants::*;
use crate::utils::vec3::*;
use cubecl::prelude::*;
use crate::simulation::accessgpuspins::get_spin;
use crate::simulation::field::{get_anisotropy_field, get_interaction_field, get_spin_orbit_field};
use crate::simulation::interactions::GPUInteractions;

#[cube]
fn llg_derivative<F: Float>(
    s: &Vec3<F>,
    b_eff: &Vec3<F>,
    alpha: F
) -> Vec3<F> {
    let prefactor:F = F::cast_from(GYROMAGNETIC_RATIO)/(F::cast_from(1.0_f32)+alpha*alpha);
    let s_cross_b = vec3_cross(s, b_eff);
    let llg_derivative = vec3_scalar_multiply(
        -prefactor, &vec3_add(
            &s_cross_b,
            &vec3_scalar_multiply(
                alpha,
                &vec3_cross(s, &s_cross_b)
            )
        )
    );
    llg_derivative
}

#[cube]
pub fn llg_rk4<F: Float>(
    interactions_data: &GPUInteractions, spins: &Array<F>,
    layers: u32, x_cells: u32, y_cells: u32,
    timestep:f32, anisotropy_strength:f32, alpha:f32,
    layer:u32, uc_x:u32, uc_y:u32, site:u32,
    num_interactions:u32,
    site_magnetic_moment:f32,
    b_ext:&Vec3<F>,
    sot_b_dl:f32, sot_b_fl:f32,
    sot_sigma:&Vec3<F>,
) -> Vec3<F> {
    let dt=F::cast_from(timestep);
    let alpha_f=F::cast_from(alpha);
    // independent of current site spin
    let interaction_field = get_interaction_field(
        interactions_data, &spins,
        layer, uc_x, uc_y, site,
        layers, x_cells, y_cells,
        site_magnetic_moment, num_interactions
    );
    let constant_field = vec3_add(&interaction_field, b_ext);

    let spin = get_spin(
        &spins,
        layer, uc_x, uc_y, site, x_cells, y_cells
    );

    let anisotropy_field_1 = get_anisotropy_field(
        &spin,
        F::cast_from(anisotropy_strength),
        site_magnetic_moment,
    );
    let sot_field_1 = get_spin_orbit_field(
        &spin, sot_sigma,
        sot_b_dl, sot_b_fl,
    );
    let varying_field_1 = vec3_add(&anisotropy_field_1, &sot_field_1);

    let b_eff_1 = vec3_add(&constant_field, &varying_field_1);
    let derivative_1 = llg_derivative(&spin, &b_eff_1, alpha_f);

    let spin_2 = vec3_add(
        &vec3_scalar_multiply(dt/F::cast_from(2.0f32), &derivative_1),
        &spin
    );
    // only need to update anisotropy field, interactions field is independent of the spin
    let anisotropy_field_2 = get_anisotropy_field(
        &spin_2,
        F::cast_from(anisotropy_strength),
        site_magnetic_moment,
    );
    let sot_field_2 = get_spin_orbit_field(
        &spin_2, sot_sigma,
        sot_b_dl, sot_b_fl,
    );
    let varying_field_2 = vec3_add(&anisotropy_field_2, &sot_field_2);

    let b_eff_2 = vec3_add(&constant_field, &varying_field_2);
    let derivative_2 = llg_derivative(&spin_2, &b_eff_2, alpha_f);

    let spin_3 = vec3_add(
        &vec3_scalar_multiply(dt/F::cast_from(2.0f32), &derivative_2),
        &spin,
    );
    let anisotropy_field_3 = get_anisotropy_field(
        &spin_3,
        F::cast_from(anisotropy_strength),
        site_magnetic_moment,
    );
    let sot_field_3 = get_spin_orbit_field(
        &spin_3, sot_sigma,
        sot_b_dl, sot_b_fl,
    );
    let varying_field_3 = vec3_add(&anisotropy_field_3, &sot_field_3);
    let b_eff_3 = vec3_add(&constant_field, &varying_field_3);
    let derivative_3 = llg_derivative(&spin_3, &b_eff_3, alpha_f);

    let spin_4 = vec3_add(
        &vec3_scalar_multiply(dt, &derivative_3),
        &spin
    );
    let anisotropy_field_4 = get_anisotropy_field(
        &spin_4,
        F::cast_from(anisotropy_strength),
        site_magnetic_moment,
    );
    let sot_field_4 = get_spin_orbit_field(
        &spin_4, sot_sigma,
        sot_b_dl, sot_b_fl,
    );
    let varying_field_4 = vec3_add(&anisotropy_field_4, &sot_field_4);
    let b_eff_4 = vec3_add(&constant_field, &varying_field_4);
    let derivative_4 = llg_derivative(&spin_4, &b_eff_4, alpha_f);

    let rk4_derivative = vec3_add(
        &vec3_add(
            &vec3_scalar_multiply(
                F::cast_from((1./6.) as f32),
                &derivative_1
            ),
            &vec3_scalar_multiply(
                F::cast_from((2./6.) as f32),
                &derivative_2
            )
        ),
        &vec3_add(
            &vec3_scalar_multiply(
                F::cast_from((2./6.) as f32),
                &derivative_3
            ),
            &vec3_scalar_multiply(
                F::cast_from((1./6.) as f32),
                &derivative_4
            )
        ),
    );

    let combined_delta = vec3_scalar_multiply(
        dt, &rk4_derivative);
    let updated_spin = vec3_norm(
        &vec3_add(&spin, &combined_delta)
    );
    updated_spin
}