use crate::simulation::unitcelldata::{UNIT_CELL_DIMENSIONS, UNIT_CELL_LAYERS_COORDINATES};
use rand::rngs::StdRng;
use rand::{Rng, RngExt};
use std::f64::consts::TAU;

fn random_unit_vector(rng: &mut impl Rng) -> [f32; 3] {
    let z: f64 = rng.random_range(-1.0..=1.0);
    let phi: f64 = rng.random_range(0.0..TAU);

    let r = (1.0 - z * z).sqrt();

    [(r * phi.cos()) as f32, (r * phi.sin()) as f32, z as f32]
}

pub fn dimensions_to_unit_cells(dims: [f32; 3]) -> [usize; 3] {
    let num_unit_cells:[f32;3] = [
        dims[0]/ UNIT_CELL_DIMENSIONS[0],
        dims[1]/ UNIT_CELL_DIMENSIONS[1],
        dims[2]/ UNIT_CELL_DIMENSIONS[2]
    ];
    [num_unit_cells[0].round() as usize, num_unit_cells[1].round() as usize, num_unit_cells[2].round() as usize]
}

fn zeros_spin_array(uc_dims: [usize;3]) -> Vec<f32>{
    let unit_cells = uc_dims[0]*uc_dims[1]*uc_dims[2];
    // 24 atoms per unit cell, 3 elements per vector
    let spins:Vec<f32> = vec![0.; unit_cells *24 *3];
    spins
}

fn set_vector(spin_array: &mut Vec<f32>, index:usize, vector:[f32;3]) {
    spin_array[index*3] = vector[0];
    spin_array[index*3+1] = vector[1];
    spin_array[index*3+2] = vector[2];
}

pub fn index_to_site(
    index: usize,
    x_cells: usize,
    y_cells: usize,
) -> (usize, usize, usize, usize) {
    let cells_per_layer = x_cells * y_cells * 4;

    let layer = index / cells_per_layer;
    let remainder = index % cells_per_layer;

    let uc_x = remainder / (y_cells * 4);
    let remainder = remainder % (y_cells * 4);

    let uc_y = remainder / 4;
    let site = remainder % 4;

    (layer, uc_x, uc_y, site)
}

pub fn use_offset(
    value:usize,
    offset:isize,
    range:usize,
) -> usize {
    let r = range as isize;
    let v = value as isize + offset;
    ((v % r + r) % r) as usize
}

pub fn site_to_index(
    layer:usize,
    uc_x:usize,
    uc_y:usize,
    site:usize,
    x_cells:usize,
    y_cells:usize,
) -> usize {
    layer*x_cells*y_cells * 4
        + uc_x*y_cells*4
        + uc_y*4
        + site
}

pub fn site_to_spacial_coords(
    layer:usize,
    uc_x:usize,
    uc_y:usize,
    site:usize,
) -> [f32;3] {
    let uc_offset:[f32;3] = [
        UNIT_CELL_DIMENSIONS[0]*uc_x as f32,
        UNIT_CELL_DIMENSIONS[1]*uc_y as f32,
        UNIT_CELL_DIMENSIONS[2]*(layer/6) as f32
    ];
    let sub_offset:[f32;3] = UNIT_CELL_LAYERS_COORDINATES[layer %6][site];

    [
        uc_offset[0]+sub_offset[0],
        uc_offset[1]+sub_offset[1],
        uc_offset[2]+sub_offset[2]
    ]
}

pub fn rng_init_spin_array(uc_dims: [usize;3], rng: &mut StdRng, mask: &[u32]) -> (Vec<f32>, u32) {
    let mut spin_array = zeros_spin_array(uc_dims);

    let atoms= spin_array.len()/3;
    let mut trueatoms = 0;
    for atom in 0..atoms {
        if mask[atom]==1u32 {
            trueatoms += 1;
            set_vector(
                &mut spin_array,
                atom,
                random_unit_vector(rng)
            );
        }
    }
    (spin_array, trueatoms)
}