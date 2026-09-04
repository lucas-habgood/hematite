use std::fs;
use std::str::FromStr;
use std::sync::OnceLock;
use cubecl::num_traits::Inv;

static CONFIG: OnceLock<Config> = OnceLock::new();

pub struct Config {
    pub iterations: u32,
    pub iters_per_plot: u32,
    pub dimensions: [f32;3],
    pub cuda: bool,
    pub timestep: f32,
    pub anisotropy_strength: f32,
    pub alpha: f32,
    pub precondition_iterations: u32,
    pub precondition_timestep: f32,
    pub applied_field: [f32;3],
    pub plot_layer: usize,
    pub interfacial_dmi_energy: f32,
    pub effective_canting_dmi_energy: f32,
    pub interface_normal: [f32;3],
    pub rng_seed: u64,
    pub site_magnetic_moment:f32,
    pub num_arrows: u32,
    pub osc_freq: f32,
    pub osc_field: [f32;3],
    pub osc_sigma: f32,
    pub initial_temperature: f32,
    pub final_temperature: f32,
    pub mask_file_name: String,
    pub stats_only: bool,
    pub output_file: String,
    pub sot_b_dl: f32,
    pub sot_b_fl: f32,
    pub sot_current_direction: [f32;3],
    pub sample_output_file: String,
    pub samples: Vec<[f32;3]>,
    pub intralayer: bool,
    pub save_at_end: bool,
    pub reload_save: bool,
    pub save_file: String,
    pub load_file: String,
    pub a_cut: bool,
    pub a_cut_reorient_arrows: bool,
}

pub fn config() -> &'static Config {
    CONFIG.get_or_init(|| read_config())
}

fn parse_iter_units(unit_string: String) -> f32 {
    let factor = match unit_string.to_uppercase().as_str() {
        "K" => 1000.,
        "M" => 1000000.,
        "B" => 1000000000.,
        &_ => 1.
    };
    factor
}

fn parse_length_units(unit_string: String) -> f32 {
    let factor = match unit_string.to_lowercase().as_str() {
        "a" => 1.0,
        "nm"=> 10.0,
        "um"=> 1.0e4,
        &_ => 1.
    };
    factor
}

fn parse_time_units(unit_string: String) -> f32 {
    let factor = match unit_string.to_lowercase().as_str() {
        "s" => 1.0,
        "ms" => 1.0e-3,
        "us" => 1.0e-6,
        "ns" => 1.0e-9,
        "ps" => 1.0e-12,
        "fs" => 1.0e-15,
        "as" => 1.0e-18,
        "zs" => 1.0e-21,
        &_ => 1.
    };
    factor
}

fn parse_energy_units(unit_string: String) -> f32 {
    let factor = match unit_string.to_lowercase().as_str() {
        "j" => 1.0,
        "mj" => 1.0e-3,
        "uj" => 1.0e-6,
        "nj" => 1.0e-9,
        "pj" => 1.0e-12,
        "fj" => 1.0e-15,
        "aj" => 1.0e-18,
        "ev" => 1.602177e-19,
        "mev" => 1.602177e-22,
        "uev" => 1.602177e-25,
        "nev" => 1.602177e-28,
        &_ => 1.
    };
    factor
}

fn parse_field_units(unit_string: String) -> f32 {
    let factor = match unit_string.to_lowercase().as_str() {
        "kt" => 1.0e3,
        "t" => 1.0,
        "mt" => 1.0e-3,
        "ut" => 1.0e-6,
        "g" => 1.0e-4,
        &_ => 1.
    };
    factor
}

fn parse_freq_units(unit_string: String) -> f32 {
    let factor = match unit_string.to_lowercase().as_str() {
        "hz" => 1.0,
        "khz" => 1.0e3,
        "mhz" => 1.0e6,
        "ghz" => 1.0e9,
        "thz" => 1.0e12,
        "phz" => 1.0e15,
        &_ => 1.
    };
    factor
}

fn split_num_suffix(s: String) -> (String, String) {
    if let Some(i) = s.rfind(|c: char| c.is_ascii_digit() || c == '.') {
        let (left, right) = s.split_at(i + 1);
        (left.to_string(), right.to_string())
    } else {
        (String::new(), s)
    }
}

fn norm_vector(v: [f32;3]) -> [f32;3] {
    let mag_squared = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    let inv_sqrt = mag_squared.inv().sqrt();
    [v[0]*inv_sqrt, v[1]*inv_sqrt, v[2]*inv_sqrt]
}


fn read_config() -> Config {
    let file_contents = fs::read_to_string("hematite.cfg");
    let mut iterations:u32=10000;
    let mut iters_per_plot:u32=1000;
    let mut dimensions:[f32;3]=[500.0, 500.0, 15.0];
    let mut cuda:bool = false;
    let mut timestep:f32=5.0e-16;
    let mut anisotropy_strength:f32=-5.0e-25;
    let mut alpha:f32=1.0;
    let mut precondition_iterations:u32=0;
    let mut precondition_timestep:f32=1.0e-17;
    let mut applied_field:[f32;3]=[0.0;3];
    let mut plot_layer:usize=0;
    let mut interfacial_dmi_energy:f32=0.0;
    let mut effective_canting_dmi_energy:f32=1.1e-22;
    let mut interface_normal:[f32;3]=[0.0,0.0,1.0];
    let mut rng_seed:u64=0;
    let mut site_magnetic_moment:f32=4.3;
    let mut num_arrows:u32=5000;
    let mut osc_freq = 0.0;
    let mut osc_field:[f32;3]=[0.0;3];
    let mut osc_sigma:f32=25.0;
    let mut initial_temperature = 0.0;
    let mut final_temperature = 0.0;
    let mut mask_file_name:String=String::from_str("mask.png").unwrap();
    let mut stats_only:bool=false;
    let mut output_file:String=String::from_str("output.txt").unwrap();
    let mut sot_b_dl:f32=0.0;
    let mut sot_b_fl:f32=0.0;
    let mut sot_current_direction:[f32;3]=[1.0,0.0,0.0];
    let mut sample_output_file:String=String::from_str("samples.txt").unwrap();
    let mut samples:Vec<[f32;3]>=vec![];
    let mut intralayer:bool=false;
    let mut save_at_end=false;
    let mut reload_save=false;
    let mut save_file=String::from_str("savestate.bin").unwrap();
    let mut load_file=String::from_str("loadstate.bin").unwrap();
    let mut a_cut:bool=false;
    let mut a_cut_reorient_arrows:bool=false;

    for line_raw in file_contents.unwrap_or(String::new()).lines() {
        // remove any whitespace
        let line: String = line_raw.chars()
            .filter(|c| !c.is_whitespace())
            .collect();
        if !line.starts_with("#")
            && let Some((key, value)) = line.split_once('=') {

            match key {
                "iterations" => {
                    let (number, suffix) = split_num_suffix(value.parse().unwrap());
                    let value_float = number.parse::<f32>().unwrap();
                    let factor = parse_iter_units(suffix);
                    iterations = (value_float*factor).round() as u32;
                }
                "iters_per_plot" => {
                    let (number, suffix) = split_num_suffix(value.parse().unwrap());
                    let value_float = number.parse::<f32>().unwrap();
                    let factor = parse_iter_units(suffix);
                    iters_per_plot = (value_float*factor).round() as u32;
                }
                "dimensions" => {
                    let split_values = value.split(",").collect::<Vec<&str>>();
                    for i in 0..3 {
                        let (number, suffix) = split_num_suffix(
                            split_values[i].parse().unwrap());
                        let value_float = number.parse::<f32>().unwrap();
                        let factor = parse_length_units(suffix);
                        dimensions[i] = value_float*factor;
                    }
                }
                "cuda" => {
                    cuda = value.parse::<bool>().unwrap();
                }
                "timestep" => {
                    let (number, suffix) = split_num_suffix(value.parse().unwrap());
                    let value_float = number.parse::<f32>().unwrap();
                    let factor = parse_time_units(suffix);
                    timestep = value_float*factor;
                }
                "anisotropy_strength" => {
                    let (number, suffix) = split_num_suffix(value.parse().unwrap());
                    let value_float = number.parse::<f32>().unwrap();
                    let factor = parse_energy_units(suffix);
                    anisotropy_strength = value_float*factor;
                }
                "alpha" => {
                    alpha = value.parse::<f32>().unwrap();
                }
                "precondition_iterations" => {
                    let (number, suffix) = split_num_suffix(value.parse().unwrap());
                    let value_float = number.parse::<f32>().unwrap();
                    let factor = parse_iter_units(suffix);
                    precondition_iterations = (value_float*factor).round() as u32;
                }
                "precondition_timestep" => {
                    let (number, suffix) = split_num_suffix(value.parse().unwrap());
                    let value_float = number.parse::<f32>().unwrap();
                    let factor = parse_time_units(suffix);
                    precondition_timestep = value_float*factor;
                }
                "applied_field" => {
                    let split_values = value.split(",").collect::<Vec<&str>>();
                    for i in 0..3 {
                        let (number, suffix) = split_num_suffix(
                            split_values[i].parse().unwrap());
                        let value_float = number.parse::<f32>().unwrap();
                        let factor = parse_field_units(suffix);
                        applied_field[i] = value_float*factor;
                    }
                }
                "plot_layer" => {
                    plot_layer = value.parse::<usize>().unwrap();
                }
                "interfacial_dmi_energy" => {
                    let (number, suffix) = split_num_suffix(value.parse().unwrap());
                    let value_float = number.parse::<f32>().unwrap();
                    let factor = parse_energy_units(suffix);
                    interfacial_dmi_energy = value_float*factor;
                }
                "effective_canting_dmi_energy" => {
                    let (number, suffix) = split_num_suffix(value.parse().unwrap());
                    let value_float = number.parse::<f32>().unwrap();
                    let factor = parse_energy_units(suffix);
                    effective_canting_dmi_energy = value_float*factor;
                }
                "interface_normal" => {
                    let split_values = value.split(",").collect::<Vec<&str>>();
                    let mut temp_array:[f32;3]=[0.0;3]; // for normalising
                    for i in 0..3 {
                        temp_array[i] = split_values[i].parse::<f32>().unwrap();
                    }
                    interface_normal = norm_vector(temp_array);
                }
                "rng_seed" => {
                    rng_seed = value.parse::<u64>().unwrap();
                }
                "site_magnetic_moment" => {
                    site_magnetic_moment = value.parse::<f32>().unwrap();
                }
                "num_arrows" => {
                    num_arrows = value.parse::<u32>().unwrap();
                }
                "oscillation_frequency" => {
                    let (number, suffix) = split_num_suffix(value.parse().unwrap());
                    let value_float = number.parse::<f32>().unwrap();
                    let factor = parse_freq_units(suffix);
                    osc_freq = value_float*factor;
                }
                "initial_temperature" => {
                    initial_temperature = value.parse::<f32>().unwrap();
                }
                "final_temperature" => {
                    final_temperature = value.parse::<f32>().unwrap();
                }
                "oscillation_field" => {
                    let split_values = value.split(",").collect::<Vec<&str>>();
                    for i in 0..3 {
                        let (number, suffix) = split_num_suffix(
                            split_values[i].parse().unwrap());
                        let value_float = number.parse::<f32>().unwrap();
                        let factor = parse_field_units(suffix);
                        osc_field[i] = value_float*factor;
                    }
                }
                "oscillation_size" => {
                    let (number, suffix) = split_num_suffix(value.parse().unwrap());
                    let value_float = number.parse::<f32>().unwrap();
                    let factor = parse_length_units(suffix);
                    osc_sigma = value_float*factor;
                }
                "mask_file_name" => {
                    mask_file_name = value.parse::<String>().unwrap();
                }
                "stats_only" => {
                    stats_only = value.parse::<bool>().unwrap();
                }
                "output_file" => {
                    output_file = value.parse::<String>().unwrap();
                }
                "sot_damping" => {
                    let (number, suffix) = split_num_suffix(value.parse().unwrap());
                    let value_float = number.parse::<f32>().unwrap();
                    let factor = parse_field_units(suffix);
                    sot_b_dl = value_float*factor;
                }
                "sot_field" => {
                    let (number, suffix) = split_num_suffix(value.parse().unwrap());
                    let value_float = number.parse::<f32>().unwrap();
                    let factor = parse_field_units(suffix);
                    sot_b_fl = value_float*factor;
                }
                "sot_current_direction" => {
                    let split_values = value.split(",").collect::<Vec<&str>>();
                    let mut temp_array:[f32;3]=[0.0;3]; // for normalising
                    for i in 0..3 {
                        temp_array[i] = split_values[i].parse::<f32>().unwrap();
                    }
                    sot_current_direction = norm_vector(temp_array);
                }
                "sample_output_file" => {
                    sample_output_file = value.parse::<String>().unwrap();
                }
                "sample" => {
                    let mut sample = [0.0;3];
                    let split_values = value.split(",").collect::<Vec<&str>>();
                    for i in 0..3 {
                        let value_float = split_values[i].parse::<f32>().unwrap();
                        sample[i] = value_float;
                    }
                    samples.push(sample);
                }
                "intralayer_interactions" => {
                    intralayer = value.parse::<bool>().unwrap();
                }
                "save_at_end" => {
                    save_at_end = value.parse::<bool>().unwrap();
                }
                "reload_save" => {
                    reload_save = value.parse::<bool>().unwrap();
                }
                "save_file" => {
                    save_file = value.parse::<String>().unwrap();
                }
                "load_file" => {
                    load_file = value.parse::<String>().unwrap();
                }
                "a_cut" => {
                    a_cut = value.parse::<bool>().unwrap();
                }
                "a_cut_reorient_arrows" => {
                    a_cut_reorient_arrows = value.parse::<bool>().unwrap();
                }
                &_ => {}
            }
        }
    }

    Config {
        iterations,
        iters_per_plot,
        dimensions,
        cuda,
        timestep,
        anisotropy_strength,
        alpha,
        precondition_iterations,
        precondition_timestep,
        applied_field,
        plot_layer,
        interfacial_dmi_energy,
        effective_canting_dmi_energy,
        interface_normal,
        rng_seed,
        site_magnetic_moment,
        num_arrows,
        osc_freq,
        osc_field,
        osc_sigma,
        initial_temperature,
        final_temperature,
        mask_file_name,
        stats_only,
        output_file,
        sot_b_dl,
        sot_b_fl,
        sot_current_direction,
        sample_output_file,
        samples,
        intralayer,
        save_at_end,
        reload_save,
        save_file,
        load_file,
        a_cut,
        a_cut_reorient_arrows,
    }
}