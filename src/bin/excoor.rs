extern crate vasp_rust;
extern crate regex;

use std::fs::File;
use std::io::BufWriter;
use std::io::prelude::*;


use regex::Regex;

fn main(){
    //cmd = "sed -n '7p' POSCAR"
    let cmd = vasp_rust::exec_cmd("sed", &["-n", "7p", "POSCAR"]);
    let total_atoms: usize = cmd.trim().split_whitespace().map(|num| num.parse::<usize>().unwrap()).sum::<usize>();


    //grep "POSITION" -A %d OUTCAR
    let total_atoms_s = format!("{}", total_atoms + 15);
    let cmd = vasp_rust::exec_cmd("grep", &["POSITION", "-A", &total_atoms_s,"OUTCAR"]);
    
    let content = cmd.split('\n').collect::<Vec<&str>>();

    let f = File::create("OUTCAR.pos").expect("Unable to write file");
    let mut output_file = BufWriter::new(f);

    let pattern = Regex::new(r"POSITION").unwrap();

    let mut idx: usize = 0;
    let mut step: usize = 0;
    while idx < content.len(){
        if !pattern.is_match(&content[idx]){
            idx += 1;
            continue;
        }
        step += 1;
        output_file.write(format!("Step: {}\n", step).as_bytes());
        output_file.write(format!("{}\n", content[idx + 1]).as_bytes());
        idx += 2;

        let end = idx + total_atoms;
        for i in idx..end{
            let line = content[i].trim().split_whitespace().map(|s| s.parse::<f64>().unwrap()).collect::<Vec<f64>>();
            let tmp = (line[3] * line[3] + line[4] * line[4] + line[5] * line[5]).sqrt();
            output_file.write(format!("{:>10.5}  {:>10.5}  {:>10.5}      {:>10.5}\n", line[0], line[1], line[2], tmp).as_bytes());
        }
        idx = end;
        output_file.write(format!("{}\n", content[idx]).as_bytes());
        output_file.write(format!("{}\n", content[idx + 10]).as_bytes());
        output_file.write(format!("{}\n\n", content[idx + 12]).as_bytes());
        idx += 14;
    }
    println!("\n    --------------------Done--------------------\n")
}