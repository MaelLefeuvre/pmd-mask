use log::LevelFilter;
use log::Level;
use env_logger::{Builder, Env, fmt::Color};
use std::io::Write;

use chrono::Local;

/// A basic static logger. This is just a quick wrapper to [`env_logger`], with additional methods.
pub struct Logger;

impl Logger {
    /// Initialize a logger with a set verbosity. 
    /// 
    /// Note that [`env_logger`] uses a static, singleton pattern, andonly a single [`Logger`] may be 
    /// initialized within a program.
    /// 
    /// # Behavior
    /// 
    /// Only four levels can be obtained from `verbosity`, and are directly tied to the external [`log::LevelFilter`] enum.
    /// A quick `verbosity` -> [`LevelFilter`](`log::LevelFilter`) table can be found below.
    /// 
    /// See the documentations of [`log::LevelFilter`] and [`log::Level`] for a more detailled explanation of each level's
    /// behavior.
    /// 
    /// | `verbosity`    | [`LevelFilter`](`log::LevelFilter`) variant |
    /// | -------------- | ------------------------------------------- |
    /// | `0`            | [`Error`](`log::LevelFilter::Error`)        |
    /// | `1`            | [`Warn`](`log::LevelFilter::Warn`)          |
    /// | `2`            | [`Info`](`log::LevelFilter::Info`)          |
    /// | `3`            | [`Debug`](`log::LevelFilter::Debug`)        |
    /// | `4..u8::MAX`   | [`Trace`](`log::LevelFilter::Trace`)        |
    /// 
    pub fn init(verbosity: u8) {
        let log_level = Self::u8_to_loglevel(verbosity);
        let env = Env::default()
            .filter("PMDMASK_LOG");

        Builder::new().filter_level(log_level)
            .format(|buf, record| {
                
                let set_intensity = record.level() == LevelFilter::Error;

                let mut arg_style = buf.style();
                arg_style.set_intense(set_intensity);


                let mut level_style = buf.style();
                let color = match record.level() {
                    Level::Error => Color::Red,
                    Level::Warn  => Color::Yellow, 
                    Level::Info  => Color::Green,
                    Level::Debug => Color::Blue,
                    Level::Trace => Color::Cyan
                };
                level_style.set_color(color).set_bold(true);

                writeln!(buf, "[{} {: <5} {}] {}",
                    Local::now().format("%Y-%m-%dT%H:%M:%S"),
                    level_style.value(record.level()),
                    record.target(),
                    arg_style.value(record.args())
                )
            })
            .parse_env(env)
            .init();
    }

    /// Convert a [`u8`] primitive into a [`log::LevelFilter`] enum variant.
    fn u8_to_loglevel(verbosity: u8) -> LevelFilter {
        match verbosity {
            0            => LevelFilter::Error,
            1            => LevelFilter::Warn,
            2            => LevelFilter::Info,
            3            => LevelFilter::Debug,
            4..= u8::MAX => LevelFilter::Trace
        }
    }

    /// Manually Set the level of a [`Logger`]
    #[cfg(test)]
    pub fn set_level(verbosity: u8) {
        log::set_max_level(Self::u8_to_loglevel(verbosity));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_level(){
        for level in 0..u8::MAX {
            Logger::set_level(level);

            let expected_level = match level {
                0           => LevelFilter::Error,
                1           => LevelFilter::Warn,
                2           => LevelFilter::Info,
                3           => LevelFilter::Debug,
                4..=u8::MAX => LevelFilter::Trace
            };

            assert_eq!(log::max_level(), expected_level);
        }
    }
}