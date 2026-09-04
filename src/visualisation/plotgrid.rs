use crate::visualisation::interpolategrid::{Point, get_bounds, make_grid};
use plotters::prelude::*;
use plotters_arrows::ThinArrow;
use crate::utils::config::config;

pub fn plot_points(
    points: &Vec<Point>,
    step_size: f32,
    filename: &str,
) {
    let (x_min, x_max, y_min, y_max) = get_bounds(points);
    let grid = make_grid(&points,step_size);

    let root = BitMapBackend::new(
        filename,
        (1000,1000)
    ).into_drawing_area();
    root.fill(&WHITE).unwrap();
    let mut chart = ChartBuilder::on(&root)
        .margin(5)
        .x_label_area_size(50)
        .y_label_area_size(50)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)
        .unwrap();
    chart.configure_mesh().x_desc("X (Angstroms)").y_desc("Y (Angstroms)").draw().unwrap();


    let dx = x_max - x_min;
    let x_steps = (dx / step_size).floor() as usize;
    let dy = y_max - y_min;
    let y_steps = (dy / step_size).floor() as usize;
    let mut cells = Vec::with_capacity((x_steps - 1) * (y_steps - 1));
    for iy in 0..y_steps - 1 {
        for ix in 0..x_steps - 1 {
            let x0 = x_min + ix as f32 * step_size;
            let x1 = x_min + (ix + 1) as f32 * step_size;
            let y0 = y_min + iy as f32 * step_size;
            let y1 = y_min + (iy + 1) as f32 * step_size;

            let value = (
                grid[ix][iy]
                    + grid[ix+1][iy]
                    + grid[ix][iy+1]
                    + grid[ix+1][iy+1]
            ) / 4.0;

            cells.push(Rectangle::new(
                [(x0, y0), (x1, y1)],
                value_to_color(value).filled(),
            ));
        }
    }

    let num_arrows:u32 = config().num_arrows;
    let arrows_per_axis = (num_arrows as f32).sqrt(); // assumes dx~=dy
    let arrow_every = points.len() as f32/ num_arrows as f32;
    let mut draw_arrows = vec![];
    for arrow in 0..num_arrows {
        let point = &points[(arrow as f32*arrow_every) as usize];
        if point.x<(x_steps-1) as f32 *step_size && point.y<(y_steps-1) as f32 *step_size {
            draw_arrows.push(ThinArrow::new(
                (point.x, point.y), (point.x+point.vx*dx/arrows_per_axis, point.y+point.vy*dy/arrows_per_axis),&BLACK
            ))
        }
    }

    chart.draw_series(cells).unwrap();
    chart.draw_series(draw_arrows).unwrap();
    root.present().unwrap();
}


fn value_to_color(value: f32) -> RGBColor {
    let (r,g,b) = coolwarm(value);


    RGBColor((r*255.0) as u8, (g*255.0) as u8, (b*255.0) as u8)
}


fn coolwarm(x: f32) -> (f32, f32, f32) {
    // Clamp x to [-1, 1]
    let x = x.clamp(-1.0, 1.0);

    let (x0, c0, x1, c1) =
    if x <= -0.3 {
        (
            -1.0,
            [0.0, 0.0, 0.6],
            -0.3,
            [0.0, 0.0, 1.0],
        )
    }
    else if x <= 0.0 {
        (
            -0.3,
            [0.0, 0.0, 1.0],
            0.0,
            [0.86, 0.86, 0.86],
        )
    }
    else if x<=0.3 {
        (
            0.0,
            [0.86, 0.86, 0.86],
            0.3,
            [1.0, 0.0, 0.0],
        )
    }
    else {
        (
            0.3,
            [1.0, 0.0, 0.0],
            1.0,
            [0.6, 0.0, 0.0],
        )
    };

    // Linear interpolation
    let t = (x - x0) / (x1 - x0);

    (
        c0[0] + t * (c1[0] - c0[0]),
        c0[1] + t * (c1[1] - c0[1]),
        c0[2] + t * (c1[2] - c0[2]),
    )
}