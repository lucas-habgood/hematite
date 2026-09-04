#[derive(Copy, Clone)]
pub struct Point {
    pub x: f32,
    pub y: f32,
    pub vx: f32,
    pub vy: f32,
    pub vz: f32,
}

fn idw_interpolate(
    x: f32,
    y: f32,
    step_size: f32,
    points: &[Point],
) -> f32 {
    let mut weighted_sum = 0.0;
    let mut weight_sum = 0.0;

    for p in points {
        let dx = x - p.x;
        let dy = y - p.y;
        let dist_sq = dx*dx + dy*dy;

        if dist_sq < 1e-10 {
            return p.vz;
        }
        if dist_sq > 10.*step_size * step_size {
            continue;
        }

        let weight = 1.0 / dist_sq.powf(0.5);

        weighted_sum += weight * p.vz;
        weight_sum += weight;
    }

    weighted_sum / weight_sum
}

fn bucket_points(
    points: &[Point],
    buckets_x: usize,
    buckets_y: usize,
) -> Vec<Vec<Vec<Point>>> {
    let (x_min, x_max, y_min, y_max) = get_bounds(&points);
    // slightly larger to not store single highest value in separate bucket
    let bucket_size_x = 1.001*(x_max - x_min)/buckets_x as f32;
    let bucket_size_y = 1.001*(y_max - y_min)/buckets_y as f32;

    let mut buckets = vec![];
    for _x in 0..buckets_x+1 {
        let mut buckets_column:Vec<Vec<Point>> = vec![];
        for _y in 0..buckets_y+1 {
            buckets_column.push(vec![]);
        }
        buckets.push(buckets_column);
    }

    for point in points {
        let bucket_x = (point.x/bucket_size_x).floor() as usize;
        let bucket_y = (point.y/bucket_size_y).floor() as usize;
        buckets[bucket_x][bucket_y].push(point.clone());
    }
    buckets
}

pub(crate) fn get_bounds(
    points: &[Point]
) -> (f32,f32,f32,f32) {
    let x_min = points.iter().map(|p| p.x)
        .fold(f32::INFINITY, f32::min);

    let x_max = points.iter().map(|p| p.x)
        .fold(f32::NEG_INFINITY, f32::max);

    let y_min = points.iter().map(|p| p.y)
        .fold(f32::INFINITY, f32::min);

    let y_max = points.iter().map(|p| p.y)
        .fold(f32::NEG_INFINITY, f32::max);
    (x_min, x_max, y_min, y_max)
}

fn near_bucket_points(
    buckets: &[Vec<Vec<Point>>],
    x: f32,
    y: f32,
    buckets_x: usize,
    buckets_y: usize,
    bucket_size_x: f32,
    bucket_size_y: f32,
) -> Vec<Point> {
    let mut near_points = vec![];

    let bucket_near_x = (x/bucket_size_x).floor() as i32;
    let bucket_near_y = (y/bucket_size_y).floor() as i32;

    for b_x_offset in -1..2 {
        for b_y_offset in -1..2 {
            let bucket_x = bucket_near_x+b_x_offset;
            let bucket_y = bucket_near_y+b_y_offset;
            if bucket_x < buckets_x as i32 && bucket_x >= 0
            && bucket_y < buckets_y as i32 && bucket_y >= 0{
                for point in buckets[bucket_x as usize][bucket_y as usize].clone() {
                    near_points.push(point)
                }
            }
        }
    }

    near_points
}

pub fn make_grid(
    points: &[Point],
    step_size: f32,
) -> Vec<Vec<f32>> {
    let (x_min, x_max, y_min, y_max) = get_bounds(points);
    let dx = x_max - x_min;
    let x_steps = (dx / step_size).floor() as usize;
    let dy = y_max - y_min;
    let y_steps = (dy / step_size).floor() as usize;

    let buckets_x = (dx/20.0).ceil() as usize;
    let buckets_y = (dy/20.0).ceil() as usize;
    let bucket_size_x = 1.001*(x_max - x_min)/buckets_x as f32;
    let bucket_size_y = 1.001*(y_max - y_min)/buckets_y as f32;
    let buckets = bucket_points(points, buckets_x, buckets_y);

    let mut grid = vec![vec![0.0; y_steps]; x_steps];

    for i in 0..x_steps {
        let x = x_min + (i as f32 * step_size);
        for j in 0..y_steps {
            let y = y_min + (j as f32 * step_size);
            let nearby_points = near_bucket_points(
                &buckets,
                x, y, buckets_x, buckets_y,
                bucket_size_x, bucket_size_y,
            );
            let interpolated = idw_interpolate(x, y, step_size, &nearby_points);
            grid[i][j] = interpolated;
        }
    }
    grid
}

