

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
