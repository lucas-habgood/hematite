use crate::utils::cpuspinarray::{site_to_index, site_to_spacial_coords};
use crate::visualisation::interpolategrid::Point;
use crate::visualisation::plotgrid::plot_points;
use std::cmp::min;
use crate::utils::config::config;

pub fn plot(
    spins: &[f32],
    uc_dims:[usize;3],
    filename: &str,
) {
    let cut_layer_to_plot =config().plot_layer;

    let (points, min_dim) = if config().a_cut {
        let layer_spins = a_cut_layer_spins(
            &spins,
            cut_layer_to_plot,
            uc_dims,
        );
        let positions = a_cut_layer_positions(
            cut_layer_to_plot, uc_dims,
        );
        let points = make_a_cut_points(
            &layer_spins, &positions
        );
        let min_dim = min(uc_dims[1], uc_dims[2]);
        (points, min_dim)
    } else {
        let layer_spins = c_cut_layer_spins(
            &spins,
            cut_layer_to_plot,
            uc_dims[0], uc_dims[1],
        );
        let positions = c_cut_layer_positions(
            cut_layer_to_plot, uc_dims[0], uc_dims[1],
        );
        let points = make_c_cut_points(&layer_spins, &positions);
        let min_dim = min(uc_dims[0], uc_dims[1]);
        (points, min_dim)
    };

    let step_size = f32::max(min_dim as f32 * 0.02, 2.5);


    plot_points(&points, step_size, filename);
}

fn make_c_cut_points(
    layer_spins: &[[f32;3]],
    positions: &[[f32;3]],
) -> Vec<Point> {
    let mut points = vec![];
    for i in 0..layer_spins.len() {
        points.push(
            Point {
                x: positions[i][0],
                y: positions[i][1],
                vx:layer_spins[i][0],
                vy:layer_spins[i][1],
                vz:layer_spins[i][2],
            }
        )
    }
    points
}

fn c_cut_layer_positions(
    layer: usize,
    x_cells: usize,
    y_cells: usize,
) -> Vec<[f32;3]> {
    let mut positions:Vec<[f32;3]>=vec![];

    for uc_x in 0..x_cells {
        for uc_y in 0..y_cells {
            for site in 0..4 {
                positions.push(site_to_spacial_coords(layer, uc_x, uc_y, site));
            }
        }
    }
    positions
}


fn c_cut_layer_spins(
    spins: &[f32],
    layer: usize,
    x_cells: usize,
    y_cells: usize,
) -> Vec<[f32;3]>{
    let mut separated = vec![];
    let start_index = layer * x_cells*y_cells * 4*3;
    let spins_to_read = x_cells*y_cells * 4;
    for i in 0..spins_to_read {
        separated.push([
            spins[start_index + 3*i],
            spins[start_index + 3*i+1],
            spins[start_index + 3*i+2],
        ]);
    }
    separated
}

fn make_a_cut_points(
    layer_spins: &[[f32;3]],
    positions: &[[f32;3]],
) -> Vec<Point> {
    let mut points = vec![];
    let (vxi,vyi,vzi) = if config().a_cut_reorient_arrows {
        (1,2,0)
    } else {
        (0,1,2)
    };
    for i in 0..layer_spins.len() {
        points.push(
            Point {
                x: positions[i][1],
                y: positions[i][2],
                vx:layer_spins[i][vxi],
                vy:layer_spins[i][vyi],
                vz:layer_spins[i][vzi],
            }
        )
    }
    points
}

fn a_cut_layer_positions(
    a_layer: usize,
    uc_dims: [usize;3],
) -> Vec<[f32;3]> {
    let mut positions = vec![];
    // a-cut, so fix x, (a axis) and go across y,z
    for uc_y in 0..uc_dims[1] {
        for c_layer in 0..uc_dims[2] * 6 {
            for site in 0..4 {
                positions.push(site_to_spacial_coords(
                    c_layer, a_layer, uc_y, site
                ));
            }
        }
    }
    positions
}

fn a_cut_layer_spins(
    spins: &[f32],
    a_layer: usize,
    uc_dims: [usize;3],
) -> Vec<[f32;3]>{
    let mut separated = vec![];
    // a-cut, so fix x, (a axis) and go across y,z
    for uc_y in 0..uc_dims[1] {
        for c_layer in 0..uc_dims[2]*6 {
            for site in 0..4 {
                let index = site_to_index(
                    c_layer, a_layer, uc_y, site, uc_dims[0], uc_dims[1]
                );
                separated.push([
                    spins[3*index],
                    spins[3*index+1],
                    spins[3*index+2],
                ]);
            }
        }
    }
    separated
}