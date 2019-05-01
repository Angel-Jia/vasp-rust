extern crate vasp_rust;

use std::env;
use vasp_rust::Data;

fn main(){
    let args: Vec<String> = env::args().collect();
    if args.len() == 1{
        println!("\nUsage: {} vasp_file1 vasp_file2 ..", args[0]);
        println!("Please try again!\n");
    }

    println!("\n############### This script converts vasp file into gview file ###############");
    println!("             ############ CONTCAR or POSCAR -> .gjf ############\n");
    
    for argv in args[1..].iter(){
        
        let mut vasp_data = Data::new();
        println!("                             Processing {}", argv);
        
        vasp_data.read_vasp(argv);
        let new_file_name;
        if argv.ends_with(".vasp"){
            new_file_name = argv.replace(".vasp", ".gjf");
        }else{
            new_file_name = argv.clone() + ".gjf";
        }
        vasp_data.write_gjf(&new_file_name);
    }
    println!("\n                     --------------- Done ---------------\n");
}