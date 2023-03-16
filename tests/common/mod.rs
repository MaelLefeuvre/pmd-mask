
/// Canonicalize paths found within a list of string arguments, such as `["--arg", "relative/path/to/file.txt"]` 
/// 
/// # Example:
/// ```
/// // Will return an array in the form of `["--arg", "/absolute/path/to/file.txt"]`
/// canon_arg!("--arg", "relative/path/to/file.txt")
/// ```
#[macro_export]
macro_rules! canon_arg {
    (flag $flag:expr)       => {[$flag]};
    ($flag:expr, $val:expr) => {[$flag, canon_arg!(val $val)] };
    (val $val:expr)         => {
        canonicalize($val)
            .expect(&format!("Failed to canonicalize {}", $val))
            .to_str()
            .expect("Invalid path")
    };
}
