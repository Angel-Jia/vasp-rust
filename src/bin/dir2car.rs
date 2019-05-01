extern crate vasp_rust;

use std::env;
use vasp_rust::Data;

fn main(){
    let args: Vec<String> = env::args().collect();
    if args.len() == 1{
        println!("\nUsage: {} vasp_file1 vasp_file2 ..", args[0]);
        println!("Please try again!\n");
    }

    println!("\n############ This script converts direct to cartesian ############");
    println!("          ############ direct -> cartesian ############\n");
    
    for argv in args[1..].iter(){
        println!("                      processing {}", argv);
        
        let mut vasp_data = Data::new();
        vasp_data.read_vasp(argv);

        let new_file_name;
        if argv.ends_with(".vasp"){
            new_file_name = argv.replace(".vasp", "-C.vasp");
        }else{
            new_file_name = argv.clone() + "-C.vasp";
        }
        vasp_data.write_vasp(&new_file_name);
    }
    println!("\n                   ---------- Done ----------\n");
}