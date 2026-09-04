use crate::simulation::unitcelldata::UNIT_CELL_DIMENSIONS;
use crate::utils::cpuspinarray::{index_to_site, site_to_spacial_coords};
use image::ImageReader;
use crate::utils::config::config;

pub fn mask_from_image(
    uc_dims: [usize; 3],
    mask_file: &str,
) -> Result<Vec<u32>, Box<dyn std::error::Error>> {
    // mask has to be u32 for passing to the gpu, there may be a nicer way but this works
    let mut mask:Vec<u32>=vec![];
    let full_path = ["masks/",mask_file].concat();
    let img = ImageReader::open(full_path)?.decode()?.to_luma8();

    // masks for a-cut are oriented differently
    let (xi, yi) = if config().a_cut {
        (1,2)
    } else {
        (0,1)
    };

    let atoms = uc_dims[0]*uc_dims[1]*uc_dims[2]*24;
    for index in 0..atoms {
        let (layer, uc_x, uc_y, site) =
            index_to_site(index, uc_dims[0], uc_dims[1]);
        let pos = site_to_spacial_coords(layer, uc_x, uc_y, site);
        let fract_x = pos[xi]/(uc_dims[xi] as f32*UNIT_CELL_DIMENSIONS[xi]);
        let fract_y = 1.0-pos[yi]/(uc_dims[yi] as f32*UNIT_CELL_DIMENSIONS[yi]);
        let (width, height) = (img.width() as f32, img.height() as f32);
        let image_x = (fract_x*width).clamp(0.0, width-1.0).round() as u32;
        let image_y = (fract_y*height).clamp(0.0, height-1.0) as u32;

        let image_pixel = img.get_pixel(image_x, image_y)[0];
        mask.push((image_pixel>=128) as u32);
    }
    Ok(mask)
}

pub fn allow_all_mask(
    uc_dims: [usize; 3],
) -> Vec<u32> {
    println!("Couldn't find mask - allowing all");
    vec![1u32; uc_dims[0]*uc_dims[1]*uc_dims[2]*24]
}