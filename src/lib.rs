#[macro_use] extern crate itertools;
extern crate regex;

use std::io::{BufReader, BufWriter};
use std::io::prelude::*;
use std::fs::File;
use std::process::Command;

use regex::Regex;
use itertools::Itertools;


pub struct Data{
    pub lattice: f64,
    pub basis: Vec<(f64, f64, f64)>,
    pub elements: Vec<String>,
    pub num_atoms: Vec<u64>,
    pub selectiveflag: String,
    pub coordinate_type: String,
    pub coordinates: Vec<(f64, f64, f64)>,
    pub selective: Vec<(char, char, char)>,
}


impl Data{
    pub fn new() -> Self{
        Data{
            lattice: 0.0,
            basis: Vec::with_capacity(3),
            elements: vec![],
            num_atoms: vec![],
            selectiveflag: "".to_string(),
            coordinate_type: "".to_string(),
            coordinates: Vec::with_capacity(100),
            selective: Vec::with_capacity(100),
        }
    }

    pub fn read_vasp(&mut self, file_name: &str){

        let f = File::open(file_name).expect(&format!("Unable to open file: {}", file_name));
        let f = BufReader::new(f);

        let input_file = f.lines().map(|line| line.unwrap()).collect::<Vec<String>>();

        self.lattice = input_file[1].trim().parse::<f64>().unwrap();
        for i in 2..5{
            let xyz: Vec<f64> = input_file[i].trim().split_whitespace().map(|num| num.parse::<f64>().unwrap()).collect();
            //let z = iter.next().unwrap().parse::<f64>().unwrap();
            self.basis.push((xyz[0], xyz[1], xyz[2]));
        }

        self.elements = input_file[5].trim().split_whitespace().map(|ele| ele.to_string()).collect();

        self.num_atoms = input_file[6].trim().split_whitespace().map(|num| num.parse::<u64>().unwrap()).collect();

        let mut offset: usize = 7;
        if input_file[offset].starts_with("S") || input_file[offset].starts_with("s"){
            self.selectiveflag.push_str("Selective dynamics");
            offset += 1;
        }

        if input_file[offset].starts_with("C") || input_file[offset].starts_with("c"){
            self.coordinate_type.push_str("Cartesian");
        }else{
            self.coordinate_type.push_str("Direct");
        }

        offset += 1;

        let total_atoms = self.num_atoms.iter().sum::<u64>();

        for i in 0..total_atoms{
            let line: Vec<&str> = input_file[i as usize + offset].trim().split_whitespace().collect();
            self.coordinates.push((line[0].parse::<f64>().unwrap(), line[1].parse::<f64>().unwrap(), line[2].parse::<f64>().unwrap()));
            if line.len() == 6{
                let s = line[3..6].iter().map(
                    |s| match s{
                        &"T" => 'T',
                        &"F" => 'F',
                        _ => panic!(format!("wrong selective flag: {} in line {}", s, i as usize + offset + 1)),
                    }).collect::<Vec<char>>();
                self.selective.push((s[0], s[1], s[2]));
            }
        }

        // 读取坐标
        if self.coordinate_type.starts_with("D") || self.coordinate_type.starts_with("d"){
            self.coordinate_type.clear();
            self.coordinate_type.push_str("Cartesian");
            dirkar(&self.basis, &mut self.coordinates);
        }
    }

    pub fn write_vasp(&self, file_name: &str){

        let f = File::create(file_name).expect(&format!("Unable to write file: {}", file_name));
        let mut output = BufWriter::new(f);
        output.write(format!("{:.3}\n", self.elements.iter().format("  ")).as_bytes());
        output.write(format!("  {:>15.10}\n", self.lattice).as_bytes());
        for i in 0..3{
            output.write(format!("  {:>15.10}  {:>15.10}  {:>15.10}\n", self.basis[i].0, self.basis[i].1, self.basis[i].2).as_bytes());
        }
        output.write(format!("{:.3}\n", self.elements.iter().format("  ")).as_bytes());
        output.write(format!("{:.3}\n", self.num_atoms.iter().format("  ")).as_bytes());
        if self.selectiveflag != "" {output.write(format!("{}\n", self.selectiveflag).as_bytes());}
        output.write(format!("{}\n", self.coordinate_type).as_bytes());

        let mut selective = self.selective.clone();;
        if self.selectiveflag != "" && self.coordinates.len() > self.selective.len(){
            let (x_f, y_f, z_f) = selective[0];
            for _ in 0..(self.coordinates.len() - self.selective.len()){
                selective.push((x_f, y_f, z_f));
            }
        }

        for (i, coord) in self.coordinates.iter().enumerate(){
            output.write(format!("  {:>15.10}  {:>15.10}  {:>15.10}", coord.0, coord.1, coord.2).as_bytes());
            match self.selectiveflag.as_str() {
                "" => output.write(b"\n"),
                _ => output.write(format!(" {} {} {}\n", selective[i].0, selective[i].1, selective[i].2).as_bytes()),
            };
        }
    }

    pub fn read_gjf(&mut self, file_name: &str){
        let f = File::open(file_name).expect(&format!("Unable to open file: {}", file_name));
        let input_file = BufReader::new(f).lines().map(|s| s.unwrap()).collect::<Vec<String>>();

        let pattern = Regex::new(r"[0-9]+\s+[0-9]+").unwrap(); //spin multiplicity in .gjf
        
        let mut space_num = 0;
        let mut idx = 0;

        loop {
            if input_file[idx].trim().len() == 0{
                space_num += 1;
                continue;
            }

            if space_num == 2 && pattern.is_match(&input_file[idx]){
                idx += 1;
                break;
            }
        }
        
        self.elements.push("".to_string());
        for i in idx..input_file.len(){
            if input_file[i].trim().len() == 0{break;}
            // element  x  y  z
            let line = input_file[i].trim().split_whitespace().collect::<Vec<&str>>();

            if line[0] == self.elements.last().unwrap(){
                (*self.num_atoms.last_mut().unwrap()) += 1;
            }else{
                self.elements.push(line[0].to_string());
                self.num_atoms.push(1);
            }

            self.coordinates.push((line[1].parse::<f64>().unwrap(), line[2].parse::<f64>().unwrap(), line[3].parse::<f64>().unwrap()));
        }
        self.elements.remove(0);
    }

    pub fn write_gjf(&self, file_name: &str){
        let f = File::create(file_name).expect("Unable to write file");
        let mut output = BufWriter::new(f);
        output.write(b"# opt freq b3lyp/6-31g\n\n");
        output.write(b"creat from vasp file\n\n");
        output.write(b"0 1\n");
        let mut total_index = 0;
        for (i, ele) in self.elements.iter().enumerate(){
            for j in 0..self.num_atoms[i]{
                output.write(format!("{:>2}    {:>15.10}  {:>15.10}  {:>15.10}\n", ele,
                        self.coordinates[total_index].0, self.coordinates[total_index].1, self.coordinates[total_index].2).as_bytes());
                total_index += 1;
            }
        }
        output.write(b"\n");
    }

    pub fn read_xyz(&mut self, file_name: &str){
        let f = File::open(file_name).expect("Unable to open file");
        let input_file = BufReader::new(f).lines().map(|s| s.unwrap().trim().to_string()).collect::<Vec<String>>();

        let total_num = input_file[0].parse::<usize>().unwrap();

        self.elements.push("".to_string());
        for i in 2..total_num{
            let line = input_file[i].split_whitespace().collect::<Vec<&str>>();
            if line[0] == self.elements.last().unwrap(){
                (*self.num_atoms.last_mut().unwrap()) += 1;
            }else{
                self.elements.push(line[0].to_string());
                self.num_atoms.push(1);
            }

            self.coordinates.push((line[1].parse::<f64>().unwrap(), line[2].parse::<f64>().unwrap(), line[3].parse::<f64>().unwrap()));
        }
        self.elements.remove(0);
    }

    pub fn write_xyz(&self, file_name: &str){
        let f = File::create(file_name).expect("Unable to write file");
        let mut output = BufWriter::new(f);
        let total_num = self.num_atoms.iter().sum::<u64>() as usize;
        let mut total_index: usize = 0;

        while total_index + total_num <= self.coordinates.len(){
            output.write(format!("{}\n", total_num).as_bytes());
            output.write(b"creat from rust\n");
            for (i, ele) in self.elements.iter().enumerate(){
                for j in 0..self.num_atoms[i]{
                    output.write(format!("{:>2}    {:>15.10}  {:>15.10}  {:>15.10}\n", ele,
                            self.coordinates[total_index].0, self.coordinates[total_index].1, self.coordinates[total_index].2).as_bytes());
                    total_index += 1;
                }
            }
        }
        output.write(b"\n");
    }
}










fn dirkar(basis: &Vec<(f64, f64, f64)>, coordinates: &mut Vec<(f64, f64, f64)>){
    for i in 0..coordinates.len(){
        let v1 = coordinates[i].0 * basis[0].0 + coordinates[i].1 * basis[1].0 + coordinates[i].2 * basis[2].0;
        let v2 = coordinates[i].0 * basis[0].1 + coordinates[i].1 * basis[1].1 + coordinates[i].2 * basis[2].1;
        let v3 = coordinates[i].0 * basis[0].2 + coordinates[i].1 * basis[1].2 + coordinates[i].2 * basis[2].2;
        coordinates[i].0 = v1;
        coordinates[i].1 = v2;
        coordinates[i].2 = v3;
    }
}

fn kardir(basis: &Vec<(f64, f64, f64)>, coordinates: &mut Vec<(f64, f64, f64)>){
    let inverse: Vec<(f64, f64, f64)> 
        = vec![(basis[1].1*basis[2].2-basis[2].1*basis[1].2, basis[2].1*basis[0].2-basis[0].1*basis[2].2, basis[0].1*basis[1].2-basis[1].1*basis[0].2),
               (basis[2].0*basis[1].2-basis[1].0*basis[2].2, basis[0].0*basis[2].2-basis[2].0*basis[0].2, basis[1].0*basis[0].2-basis[0].0*basis[1].2),
               (basis[1].0*basis[2].1-basis[2].0*basis[1].1, basis[2].0*basis[0].1-basis[0].0*basis[2].1, basis[0].0*basis[1].1-basis[1].0*basis[0].1)];
    let omega: f64 = basis[0].0*basis[1].1*basis[2].2 + basis[0].1*basis[1].2*basis[2].0 + basis[0].2*basis[1].0*basis[2].1 -
                basis[0].2*basis[1].1*basis[2].0 + basis[1].2*basis[2].1*basis[0].0 + basis[2].2*basis[0].1*basis[1].0;
    
    let inverse: Vec<(f64, f64, f64)>
        = vec![(inverse[0].0 / omega, inverse[0].1 / omega, inverse[0].2 / omega),
               (inverse[1].0 / omega, inverse[1].1 / omega, inverse[1].2 / omega),
               (inverse[2].0 / omega, inverse[2].1 / omega, inverse[2].2 / omega)];

    for i in 0..coordinates.len(){
        let v1: f64 = coordinates[i].0 * basis[0].0 + coordinates[i].1 * basis[1].0 + coordinates[i].2 * basis[2].0;
        let v2: f64 = coordinates[i].0 * basis[0].1 + coordinates[i].1 * basis[1].1 + coordinates[i].2 * basis[2].1;
        let v3: f64 = coordinates[i].0 * basis[0].2 + coordinates[i].1 * basis[1].2 + coordinates[i].2 * basis[2].2;

        // move atoms to primative cell
        coordinates[i].0 = v1 + 100 as f64 - ((v1 + 100 as f64) as i64) as f64;
        coordinates[i].1 = v2 + 100 as f64 - ((v2 + 100 as f64) as i64) as f64;
        coordinates[i].2 = v3 + 100 as f64 - ((v3 + 100 as f64) as i64) as f64;
    }
}

pub fn exec_cmd(cmd: &str, args: &[&str]) -> String{
    let output = Command::new(cmd).args(args).output().unwrap();
    if !output.status.success(){
        println!("error with execute cmd: {} {}", cmd, args.iter().format("  "));
        panic!("error message: {}", String::from_utf8_lossy(&output.stderr));
    }else{
        String::from_utf8(output.stdout).unwrap()
    }
}


// fn main(){
//     let mut lattice: f64 = 0.0;
//     let mut basis: Vec<(f64, f64, f64)> = Vec::new();
//     let mut elements: Vec<String> = Vec::new();
//     let mut num_atoms: Vec<u64> = Vec::new();
//     let mut selectiveflag: String = String::new();
//     let mut coordinate_type: String = String::new();
//     let mut coordinates: Vec<(f64, f64, f64)> = Vec::new();
//     let mut selective: Vec<(String, String, String)> = Vec::new();
//     read_vasp("POSCAR", &mut lattice, &mut basis, &mut elements, &mut num_atoms, &mut selectiveflag, &mut coordinate_type, &mut coordinates, &mut selective);
//     write_gjf("POSCAR.gjf", &elements, &num_atoms, &coordinates);
// }