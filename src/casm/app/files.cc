#include "files.hh"

#include <cstring>

#include "casm_functions.hh"
#include "casm/misc/algorithm.hh"
#include "casm/casm_io/FileEnumerator.hh"

namespace CASM {


  // ///////////////////////////////////////
  // 'files' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int files_command(
    int argc,
    char *argv[],
    PrimClex *_primclex,
    std::ostream &sout,
    std::ostream &serr) {
    
    po::variables_map vm;
    fs::path out_path;
    std::string settings;
    std::vector<std::string> calc;
    bool gz_flag(false);

    /// Set command line options using boost program_options
    po::options_description desc("'casm files' usage");
    desc.add_options()
    ("help,h", "Write help documentation")
    
    ("settings", 
      po::value<std::string>(&settings)->default_value("curr"), 
      "Which settings to include. "
      "One of: 'curr' (default), 'all'.")
    
    ("calc", 
      po::value<std::vector<std::string> >(&calc)->multitoken(), 
      "Which calculation files to include. "
      "May be zero (default) or more of: 'settings', 'status', 'all'.")
    
    ("archive,a", 
      po::value<std::vector<std::string> >(&calc)->multitoken(), 
      "Which calculation files to include. "
      "May be zero (default) or more of: 'settings', 'status', 'all'.")
    
    ("gzip,z", 
      po::value(&gz_flag)->zero_tokens(), 
      "Write gzipped output file.")
    
    ("output,o", 
      po::value<fs::path>(&out_path), 
      "Name for output file. Use STDOUT to print results without extra messages. ");
    
    
    try {
      po::store(po::parse_command_line(argc, argv, desc), vm); // can throw
      
      if(!vm.count("help")) {
        
        std::vector<std::string> allowed;
        
        allowed = std::vector<std::string>{"", "curr", "all"};
        if(!contains(allowed, settings)) {
          std::cerr << "Error in 'casm files'. '" << settings << "' is not a valid option to --settings." << std::endl;
          std::cerr << "  Valid options are: 'curr', 'all'" << std::endl;
          return ERR_INVALID_ARG;
        }
        
        allowed = std::vector<std::string>{"settings", "status", "all"};
        for(auto it=calc.begin(); it!=calc.end(); ++it) {
          if(!contains(allowed,*it)) {
            std::cerr << "Error in 'casm files'. '" << *it << "' is not a valid option to --calc." << std::endl;
            std::cerr << "  Valid options are: 'settings', 'status', 'all'" << std::endl;
            return ERR_INVALID_ARG;
          }
        }
      }
      
      /** --help option
       */
      if(vm.count("help")) {
        std::cout << "\n";
        std::cout << desc << std::endl;


        std::cout << "DESCRIPTION  \n"
                     "    Enumerate files used by this CASM project\n"
                     "    - If --all given, include files from all settings (bset, \n"
                     "      calctype, ref, clex, eci). Otherwise only include files\n"
                     "      relevant to the current settings.\n"
                     "    - If --calc given, include training data files. This option\n"
                     "      accepts one or more of: 'settings', 'status', or 'all'.\n"  
                     "    - If --calc 'settings' given, include files from the       \n"   
                     "      'training_data' settings directories.\n"
                     "    - If --calc 'status' given, include 'properties.calc.json' \n"
                     "      and 'status.json' files.\n"
                     "    - If --calc 'all' given, include all files in the          \n"
                     "      'training_data' directory, recursively.\n\n"
                     
                     "    Examples:\n"
                     "      casm files -o files.txt\n"
                     "      - Output basic project files for the current settings.\n\n"
        
                     "      tar -czvf proj.tar.gz `casm files -o STDOUT` \n"
                     "      tar -czvf proj.tar.gz -T files.txt \n"
                     "      - Use tar to create an archive of your CASM project. \n\n"
                     
                     "      rsync -avP --relative `casm files -o STDOUT` destination -n \n"
                     "      rsync -avP --relative --files-from=files.txt destination -n \n"
                     "      - Use rsync to syncronize casm project files to 'destination' \n"
                     "      - \"-n\" option is for a \"dry-run\". Remove it when ready \n"
                     "        to do the transfer.\n";
                     
        return 0;
      }
      
      

      po::notify(vm); // throws on error, so do after help in case
      // there are any problems

    }
    catch(po::error &e) {
      std::cerr << desc << std::endl;
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(std::exception &e) {
      std::cerr << desc << std::endl;
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      return ERR_UNKNOWN;

    }
    
    // find project root
    fs::path root;
    if(!_primclex) {
      root = find_casmroot(fs::current_path());
      if(root.empty()) {
        serr << "Error in 'casm files': No casm project found." << std::endl;
        return ERR_NO_PROJ;
      }
    }
    else {
      root = _primclex->get_path();
    }
    
    auto check_gz = [ = ](fs::path p) {
      if(p.extension() == ".gz" || p.extension() == ".GZ") {
        return true;
      }
      return false;
    };
    
    if(check_gz(out_path)) {
      gz_flag = true;
    }

    // set output_stream: where the file paths are written
    std::unique_ptr<std::ostream> uniq_fout;
    std::ostream &output_stream = make_ostream_if(vm.count("output"), sout, uniq_fout, out_path, gz_flag);
    
    // set status_stream: where query settings and PrimClex initialization messages are sent
    std::ostream &status_stream = (out_path.string() == "STDOUT") ? serr : sout;

    // If '_primclex', use that, else construct PrimClex in 'uniq_primclex'
    // Then whichever exists, store reference in 'primclex'
    std::unique_ptr<PrimClex> uniq_primclex;
    PrimClex &primclex = make_primclex_if_not(_primclex, uniq_primclex, root, status_stream);
    
    std::vector<fs::path> files;
    auto result = std::back_inserter(files);
    FileEnumerator enumerate(primclex, settings=="all");
    
    result = enumerate.basic_files(result);
    result = enumerate.bset_files(result);
    result = enumerate.reference_files(result);
    result = enumerate.eci_files(result);
    
    bool calc_settings = contains(calc,"settings");
    bool calc_status = contains(calc,"status");
    bool calc_all = contains(calc,"all");
    
    if(calc_settings && !calc_all) {
      result = enumerate.calc_settings_files(result);
    }
    if(calc_status && !calc_all) {
      result = enumerate.calc_status_files(result);
    }
    if(calc_all) {
      result = enumerate.all_calc_files(result);
    }
  
    for( auto f : files) {
      output_stream << f.string() << "\n";
    }
    
    return 0;
  };

}

