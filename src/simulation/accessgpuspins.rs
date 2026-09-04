use cubecl::prelude::*;
use crate::utils::vec3::*;

#[cube]
pub fn get_array_index(
    layer:u32,
    uc_x:u32,
    uc_y:u32,
    site:u32,
    x_cells:u32,
    y_cells:u32,
) -> usize {
    (
        (layer as usize) * x_cells as usize * y_cells as usize * 4
            + (uc_x as usize)*(y_cells as usize)*4
            + (uc_y as usize)*4
            + site as usize
    )*3 // add +1, +2 for y,z components
}

#[cube]
pub fn get_spin<F: Float>(
    spins:&Array<F>,
    layer:u32,
    uc_x:u32,
    uc_y:u32,
    site:u32,
    x_cells:u32,
    y_cells:u32,
) -> Vec3<F> {
    let array_index = get_array_index(layer, uc_x, uc_y, site, x_cells, y_cells);
    Vec3::<F> {
        x: spins[array_index],
        y: spins[array_index+1],
        z: spins[array_index+2]
    }
}

#[cube]
pub fn use_offset(
    value:u32,
    offset:i32,
    range:u32,
) -> u32 {
    let r = range as i32;
    let v = value as i32 + offset;
    ((v % r + r) % r) as u32
}

#[cube]
pub fn pos_from_index(
    index:u32,
    x_cells:u32,
    y_cells:u32,
) -> (u32, u32, u32, u32) {
    let cells_per_layer = x_cells as usize * y_cells as usize * 4;

    let layer = index as usize / cells_per_layer;
    let remainder = index as usize % cells_per_layer;

    let uc_x = remainder / (y_cells as usize * 4);
    let remainder = remainder % (y_cells as usize * 4);

    let uc_y = remainder / 4;
    let site = remainder % 4;

    (layer as u32, uc_x as u32, uc_y as u32, site as u32)
}