use std::f32::consts::PI;
use crate::utils::cpuspinarray::{site_to_index, use_offset};

fn stencil_in_mask (
    stencil: [[usize; 2]; 5],
    uc_dims:[usize; 3],
    mask: &[u32],
) -> bool {
    for [ux,uy] in stencil.iter() {
        let index = site_to_index(
            0, *ux, *uy, 0,
            uc_dims[0], uc_dims[1]);
        if mask[index]==0 {
            return false
        }
    }

    true
}

fn get_stencil (
    uc_x: usize,
    uc_y: usize,
    uc_dims:[usize; 3],
) -> [[usize; 2]; 5] {
    let left_uc_x = use_offset(uc_x, -1, uc_dims[0]);
    let right_uc_x = use_offset(uc_x, 1, uc_dims[0]);
    let down_uc_y = use_offset(uc_y, -1, uc_dims[1]);
    let up_uc_y = use_offset(uc_y, 1, uc_dims[1]);

    [
        [uc_x, uc_y],
        [left_uc_x, uc_y],
        [right_uc_x, uc_y],
        [uc_x, down_uc_y],
        [uc_x, up_uc_y],
    ]
}

fn get_neel (
    spins: &[f32],
    uc_x: usize,
    uc_y: usize,
    uc_dims:[usize; 3],
) -> [f32; 3] {
    let layer_zero_index = 3*site_to_index(
        0, uc_x, uc_y, 0,
        uc_dims[0], uc_dims[1]);
    let layer_one_index = 3*site_to_index(
        1, uc_x, uc_y, 1,
        uc_dims[0], uc_dims[1]);

    let layer_zero_spin:[f32;3] = [
        spins[layer_zero_index],
        spins[layer_zero_index + 1],
        spins[layer_zero_index + 2],
    ];
    let layer_one_spin:[f32;3] = [
        spins[layer_one_index],
        spins[layer_one_index + 1],
        spins[layer_one_index + 2],
    ];

    let diff:[f32;3] = [
        layer_zero_spin[0] - layer_one_spin[0],
        layer_zero_spin[1] - layer_one_spin[1],
        layer_zero_spin[2] - layer_one_spin[2],
    ];
    let magnitude:f32 = diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
    let inv_sqrt_mag = 1.0 / magnitude.sqrt();

    [
        diff[0]*inv_sqrt_mag,
        diff[1]*inv_sqrt_mag,
        diff[2]*inv_sqrt_mag,
    ]
}

fn charge_density (
    spins: &[f32],
    stencil: [[usize; 2]; 5],
    uc_dims:[usize; 3],
) -> f32 {
    let mut neel_vectors:[[f32;3];5] = [[0.0;3];5];
    for stencil_index in 0..5 {
        let uc_x = stencil[stencil_index][0];
        let uc_y = stencil[stencil_index][1];
        neel_vectors[stencil_index] = get_neel(spins, uc_x, uc_y, uc_dims);
    }


    // need the neel vector at the site, and the x,y derivatives of it

    // don't need to consider spacing in x y direction, it cancels out with area integration
    let del_x = [
        (neel_vectors[2][0] - neel_vectors[1][0])/2.0,
        (neel_vectors[2][1] - neel_vectors[1][1])/2.0,
        (neel_vectors[2][2] - neel_vectors[1][2])/2.0,
    ]; // partial in x
    let del_y = [
        (neel_vectors[4][0] - neel_vectors[3][0])/2.0,
        (neel_vectors[4][1] - neel_vectors[3][1])/2.0,
        (neel_vectors[4][2] - neel_vectors[3][2])/2.0,
    ]; // partial in y

    // now perform L . (del_x L  x  del_y L)
    let cross_dels:[f32;3] = [
        del_x[1]*del_y[2] - del_x[2]*del_y[1],
        del_x[2]*del_y[0] - del_x[0]*del_y[2],
        del_x[0]*del_y[1] - del_x[1]*del_y[0],
    ];

    let dot = neel_vectors[0][0]*cross_dels[0]
        + neel_vectors[0][1]*cross_dels[1]
        + neel_vectors[0][2]*cross_dels[2];

    dot/(4.0*PI)
}

pub fn get_topological_charge(
    spins: &[f32],
    uc_dims:[usize;3],
    mask: &[u32],
) -> f32 {
    let mut charge:f32 = 0.0;
    // use layer 0 and 1 for neel vector grid, then
    // only compute per unit cell as these are on a fixed grid
    // layer 0 site 0 <=> layer 1 site 1 same xy position in uc
    for uc_x in 0..uc_dims[0] {
        for uc_y in 0..uc_dims[1] {
            let stencil = get_stencil(uc_x, uc_y, uc_dims);
            if stencil_in_mask(stencil, uc_dims, mask) {
                charge += charge_density(spins, stencil, uc_dims);
            }
        }
    }
    charge
}