use crate::utils::config::config;
use std::fs::{self, OpenOptions};
use std::io;

pub fn create_files_folders() -> io::Result<()> {
    // reset the plots directory
    let _ = fs::remove_dir_all("plots");
    fs::create_dir_all("plots")?;

    // do not reset savestates and masks, but ensure existence
    fs::create_dir_all("savestates")?;
    fs::create_dir_all("masks")?;

    // config file must exist, even if empty, so create but don't delete contents
    OpenOptions::new()
        .create(true)
        .append(true)
        .open("hematite.cfg")?;

    // output and samples should be erased before starting a run
    OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open(config().output_file.clone())?;
    OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open(config().sample_output_file.clone())?;

    Ok(())
}