//use flate2::read::GzDecoder; // flate2 bug on concatated files while decomression, commented!
//use rust_htslib::bgzf::Reader; // omit flate to once success full
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io;
use std::io::prelude::*;
//use fastq::Parser; // not used directly but keep here for future reference
use fastq::Record;
use std::str;
//use std::path::Path;
use std::path::PathBuf;
use clap::{arg, Command, value_parser }; // ArgAction, Command};
use regex::Regex;
//use regex::RegexBuilder;

// Capacity 
const CAPACITY: usize = 10240; // will be used on flate2::GzDecoder,commented now!
fn check_inputfiles(inputfilename:&PathBuf)->Result<PathBuf,io::Error>{
    //check if path exists
    if ! inputfilename.exists() { return Err(std::io::Error::new(io::ErrorKind::NotFound,"File not found"))  }
    // check if parent directory writable !!!
    let inputdirectory=inputfilename.parent().unwrap();
    // give this error while creating output file
    // if ! std::fs::metadata(inputfilename).unwrap().permissions().readonly() { return(Err(std::io::Error::new(io::ErrorKind::PermissionDenied,"Can't write to directory")))       }
    let infilename = inputfilename.file_name().unwrap();
    let infilestem=inputfilename.file_stem().unwrap();
    let infilestem2=PathBuf::from(infilestem).file_stem().unwrap().to_owned();
    let mut outfilename=String::new();
    outfilename.push_str(infilestem2.to_str().unwrap());
    outfilename.push_str("_renamed.fastq.gz");
    let mut outputbuffer=PathBuf::from(&inputdirectory);
    outputbuffer.push(&outfilename);
    //check if filename contains "gz"
    match infilename.to_str() {
        Some(s) => {
                        if ! s.ends_with("gz"){ return Err(std::io::Error::new(io::ErrorKind::InvalidInput,"Not Gzip file")) }
                },
        _=>{},
    }
    Ok(outputbuffer)
    
}
fn convert_fastq(inputfilename:&PathBuf , outputfilename:&PathBuf ) ->Result<(),io::Error>{
    // Input values:
    //
    // While using gzip decoder from flate2 
    let in_fh = std::fs::File::open(inputfilename).unwrap();
    let in_gz = MultiGzDecoder::new(in_fh);
    let in_buf = io::BufReader::with_capacity(CAPACITY, in_gz);

    let regexbuilder=Regex::new(r"([A-Z]\d+)L(\d)C(\d\d\d)R(\d\d\d)(\d+)\/(\d)$").unwrap();
    //Output values
    //let out_filename = "output.fastq.gz";
    let out_fh = std::fs::File::create(outputfilename)?;
    let out_gz = GzEncoder::new(out_fh, Compression::default());
    let mut out_buf = io::BufWriter::new(out_gz);
    let mut header_buffer_string=String::with_capacity(200);
    
    // Read using the fastq::Parser
    let parser = fastq::Parser::new(in_buf);
    let mut readcount=0;
    parser.each( |record| {
        readcount+=1;
        let mut c_record = record.to_owned_record();
        // Update to the header
        //let header:&str = str::from_utf8(&c_record.head).unwrap();
        header_buffer_string.clear();
        header_buffer_string.push_str(str::from_utf8(&c_record.head).unwrap());
        //match mgi_readheader2_illuminaheader(&header_buffer_string,&regexbuilder){
        match mgi_readhedaer2_illuminahederNoRegex(&header_buffer_string){
        //match mgi_readhedaer2_illuminahederNoRegex(&header_buffer_string,&lindex,cindex,rindex,pairindex ){
            Ok(new_header) =>{
                            //println!("{}",&new_header.capacity());
                            c_record.head=new_header.as_bytes().to_vec();
                            match c_record.write(&mut out_buf){
                                    Ok(_) => true,// if writes success fully continue parsing
                                    _ => false // if write not successfull  stop parsing 
            }
            
            }
            _ => false
        }


        // write to the output buffer
    }
        ).expect("Invalid FASTQ file");
    if readcount==0 {
        // remove the output file
        std::fs::remove_file(outputfilename)?;
        return Err(std::io::Error::new(io::ErrorKind::InvalidInput,"0 records parsed in fast file"))
    }

    // Flush the remaing of the buffer to the file before exit.
    out_buf.flush().expect("Can't buffer flush to the file");

    println!("Reads {} parsed in file: {:?}",readcount,inputfilename);
    Ok(())
}

fn mgi_readheader2_illuminaheader(inputstring: &str,regex_builder: &regex::Regex) -> Result<String, String> {
    // Lets try sring with capacity
    //let mut output=String::new();
    let mut output=String::with_capacity(100);
    output.push_str("M00001:1:");
    // delete from here
    let rindex=inputstring.rfind("R").unwrap();
    println!("{}",&inputstring[rindex..rindex+2]);
    println!("{:?}",rindex);
    // delete to here

    //if let Some(captures) = Regex::new(r"([A-Z]\d+)L(\d)C(\d\d\d)R(\d\d\d)(\d+)\/(\d)$").unwrap().captures(inputstring)
    if let Some(captures) = regex_builder.captures(inputstring)
    {
        //fc = captures[1].to_string();
        output.push_str(&captures[1]);
        output.push(':');


        //l = captures[2].to_string();
        output.push_str(&captures[2]);
        output.push(':');

        //tile = captures[5].to_string();
        //tile = tile.trim_start_matches('0').to_string();
        output.push_str(&captures[5].trim_start_matches('0'));
        output.push(':');

        //c = captures[3].to_string();
        //c = c.trim_start_matches('0').to_string();
        output.push_str(&captures[3].trim_start_matches('0'));
        output.push(':');


        //r = captures[4].to_string();
        //r = r.trim_start_matches('0').to_string();
        output.push_str(&captures[4].trim_start_matches('0'));

        output.push(' ');

        //pair=captures[6].to_string();
        output.push_str(&captures[6]);
        output.push_str(":N:0:1");
        //let output=format!("M00001:1:{}:{}:{}:{}:{} {}:N:0:1", fc, l, tile, c, r,pair);    Ok(output)
        Ok(output)
    }else{
        Err(format!("Read name: {}\n\t format does not match to pattern '([A-Z]\\d+)L(\\d)C(\\d\\d\\d)R(\\d\\d\\d)(\\d+)\\/(\\d)$'",inputstring))
    }

}
fn mgi_readhedaer2_illuminahederNoRegex( inputstring: &str )->Result<String,String>{
    let mut output=String::with_capacity(100);
    let lindex=inputstring.rfind("L").unwrap();
    let cindex=inputstring.rfind("C").unwrap();
    let rindex=inputstring.rfind("R").unwrap();
    let pairindex=inputstring.rfind("/").unwrap();
    /*
    println!("read ilm: {}", &inputstring);
    //println!("Capture1: {}",&inputstring[0..lindex]);
    println!("Capture2: {}", &inputstring[lindex+1..cindex]);
    println!("Capture5: {}",&inputstring[rindex+4..pairindex].trim_start_matches('0'));
    println!("Capture3: {}",&inputstring[cindex+1..rindex].trim_start_matches('0'));
    println!("Capture4: {}",&inputstring[rindex+1..rindex+4].trim_start_matches('0'));
    println!("Captures6: {}", &inputstring[pairindex+1..] );
    */
    output.push_str("M00001:1:");
    // push flowcell
    output.push_str(&inputstring[0..lindex]);
    //push lane
    output.push(':');
    output.push_str(&inputstring[lindex+1..cindex]);
    output.push(':');
    // push read number
    output.push_str(&inputstring[rindex+4..pairindex].trim_start_matches('0'));
    output.push(':');
    // push the c
    output.push_str(&inputstring[cindex+1..rindex].trim_start_matches('0'));
    output.push(':');
    // push the r
    output.push_str(&inputstring[rindex+1..rindex+4].trim_start_matches('0'));
    output.push(' ');

    // push  the pair
    output.push_str(&inputstring[pairindex+1..]);
    output.push_str(":N:0:1");

    //println!("Current ouput: {}",&output);
    Ok(output)
}

fn main() {
    // Parse cli arguments with clap
    let matches= Command::new("mgi_fastq_converter")
        .version("1.1")
        .author("Ibrahim K. <kisakesenhi@gmail.com>")
        .about("Converts the readname format to illumina readname format")
        .arg_required_else_help(true)
        .arg(
            arg!( -f --fastq <FILE> "Fastq files fastq.gz"
                )
            .takes_value(true)
            .multiple(true)
            .hide_short_help(false)
            .value_parser(value_parser!(PathBuf)),
            )
        .get_matches();
    // use with getraw and convert into iter
    if matches.contains_id("fastq"){ // another if
    if matches.value_source("fastq").expect("checked contains_id") == clap::ValueSource::CommandLine {
        let fastqfiles_itr = matches.get_raw("fastq")
            .expect("`fastq` is required")
            .into_iter();
        for fq in fastqfiles_itr{
            // Check fastq file and return proper error!
            match check_inputfiles(&PathBuf::from(fq)){
                Ok(outputbuffer) =>{
                                match convert_fastq(&PathBuf::from(fq),&outputbuffer) {
                                    Ok(_) =>{},
                                    Err(e) => {eprintln!("Error during conversion of file {:?} with error: {}",fq,e)},
                                    }
                            },
                Err(e)=> eprintln!("Failed to parse file:{:?} with error: {}",fq,e),
                }
            }   
        }
    }
}
