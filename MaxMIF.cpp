/* Author:  Yingnan Hou<houyingnan3@gmail.com> 2017
 * Usage: This is MaxMIF package. Use, redistribution, modify without limitations
 */
/*
 * Modified/enhanced by zmi, 2019
 */
/***********************************************************************/

#include "MaxMIF.h"
#include "argagg.hpp" // for cli

// DONE #include "argagg.hpp" // zmi: required for CLI
// DONE add json config file for default input files : find json library
// DONE add -o output parameter for default output folder
// DONE separate input as input matrix M or input Frequency
/***********************************************************************/
// defined CLI
// maxMIF -i MutationMatrixOrFrequencyFilePath -f true  -d true -n HumanNet1/STRING10 -o OutputFilePath
// maxMIF --input MutationMatrixOrFrequencyFilePath --frequency true --default true --network HumanNet1/STRING10 --output OutputFilePath

// maxMIF -i MutationMatrixOrFrequencyFilePath -f true -d false -n NetworkName -p NetworkFilePath -o OutputFilePath
// maxMIF --input MutationMatrixOrFrequencyFilePath --frequency true --default false --network NetworkName -path NetworkFilePath --output OutputFilePath

const string Version{"1.1,  2019 zmi"};

struct CLI {
    string sInputFile;
    bool bMutationFrequency;
    bool bDefaultNetwork;
    string sNetworkName;
    string sNetworkFile;
    string sOutputFile;
};

bool parseInputArguments(int argc, char **pString, const string &ver, CLI &cli);

template<class Element, class Container>
bool in_array(const Element &element, const Container &container) {
    return std::find(std::begin(container), std::end(container), element)
           != std::end(container);
}

int main(int argc, char **argv) {
    CLI cli;
    parseInputArguments(argc, argv, Version, cli); // get cli parameters
    MaxMIF maxmif;

    // DONE change paths and logic of main using cli struct

    cout << string(120, '-') << endl;
    cout << "MaxMIF: method for identifying cancer driver genes trough maximal mutational impact function" << endl;
    cout << "        using PPI (protein-protein interaction) networks HumanNet v1 and/or STRINGv10" << endl << endl;
    cout << "        original source code version 20190916,  extended and modified by MZ in 2019" << endl;
    cout << string(120, '-') << endl;

    // Read CancerGenes only if needed
    if (MaxMIF::bUseCancerReferenceData) {
        maxmif.ReadCancerGene();
    }
    /* Read in the reference gene sets */
    maxmif.ReadGeneID();

    /* Read in the mutation data */

    cout << endl << "Read Mutation data: " << cli.sInputFile << "   ...";
    if (!maxmif.ReadMutData(cli.sInputFile, cli.bMutationFrequency)) {
        cout << "Error during reading Mutation data. Exit..." << endl;
        return -1;
    };
    /* Compute Mutation data  */
    if(!cli.bMutationFrequency){
        cout << endl << "Compute Mutation data ..." << endl;
        maxmif.ComputeMut();
    }

    /* choose network  */
    string networkPath;
    cout << endl << "Read Network data: " << cli.sNetworkName << " ...";
    if (!maxmif.ReadNetwork(cli.sNetworkFile, cli.sNetworkName)) {
        cout << "Error during reading Network data. Exit..." << endl;
        return -1;
    }
    /* Compute MaxMIF Score */
    cout << endl << "Compute MaxMIF Score ..." << endl;
    maxmif.ComputeRank();
    /* Output MaxMIF Score */
    cout << endl << "Output file: " << cli.sOutputFile << endl;
    maxmif.OutputRank(cli.sOutputFile);

    /* Clean up */
    maxmif.NetworkClear();
    maxmif.MutDataClear();
    cout << endl << "SUCCESS" << endl;
    return 0;
}

bool parseInputArguments(int argc, char **pString, const string &ver, CLI &cli) {

    ostringstream usage;
    usage
            << endl
            << "Usage: " << pString[0]
            << " -i MutationMatrixOrFrequencyFilePath -f true  -d true -n HumanNet1/STRING10 -o OutputFilePath" << endl
            << "Usage: " << pString[0]
            << " -i MutationMatrixOrFrequencyFilePath -f true -d false -n NetworkName -p NetworkFilePath -o OutputFilePath"
            << endl
            << endl;

    argagg::parser argparser{{
                                     {"help", {"-h", "--help"},
                                             "shows this help message", 0},
                                     {"input", {"-i", "--input"},
                                             "input mutation frequency/matrix ( see -f/--frequency for type)", 1},
                                     {"frequency", {"-f", "--frequency"},
                                             " true if input file is a mutation frequency file and false if it is a mutation matrix", 1},
                                     {"default", {"-d", "--default"},
                                             " true if default predefined networks HumanNet1 or STRING10 should be used", 1},
                                     {"network", {"-n", "--network"},
                                             "name of network: either default predefined HumanNet1 / STRING10, or a new user-defined network name  together with the option -p/--path of this network", 1},
                                     {"path", {"-p", "--path"},
                                             "full path to user-defined network", 1},
                                     {"output", {"-o", "--output"},
                                             "full path to output file", 1},
                                     {"version", {"-v", "--version"},
                                             "show program version", 0},
                             }};

    argagg::parser_results args;
    try {
        args = argparser.parse(argc, pString);
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        exit(1);
    }

    if (args["help"]) {
        std::cerr << string(120, '-') << endl;
        std::cerr << pString[0] << ": " << endl;
        std::cerr << "method for identifying cancer driver genes trough maximal mutational impact function" << endl;
        std::cerr << "using PPI (protein-protein interaction) networks HumanNet v1 and/or STRINGv10" << endl << endl;
        std::cerr << "original source code version 20190916,  extended and modified by MZ in 2019" << endl;
        std::cerr << string(120, '-') << endl;
        std::cerr << argparser;
        std::cerr << string(120, '-') << endl;
        std::cerr << usage.str();
        std::cerr << string(120, '-') << endl;
        exit(0);
    }

    if (args["version"]) {
        cerr << pString[0] << " " << ver << endl;
        exit(0);
    }


    // TODO fill in cli struct;
    if (!args["input"]) {
        cerr << args.program << ": '--input' ('-i') option required" << endl;
        exit(1);
    }
    if (!args["frequency"]) {
        cerr << args.program << ": '--frequency' ('-f') option required" << endl;
        exit(1);
    }
    if (!args["default"]) {
        cerr << args.program << ": '--default' ('-d') option required" << endl;
        exit(1);
    }
    if (!args["network"]) {
        cerr << args.program << ": '--network' ('-n') option required" << endl;
        exit(1);
    }
    if (!args["output"]) {
        cerr << args.program << ": '--output' ('-o') option required" << endl;
        exit(1);
    }

    string sFalse[6] = {
            "false",
            "FALSE",
            "False",
            "f",
            "F",
            "0"};

    string sTrue[6] = {
            "true",
            "TRUE",
            "True",
            "t",
            "T",
            "1"};

    string sDefaultNetworks[2] = {
            "STRING10",
            "HumanNet1"};

    if (in_array(args["default"].as<std::string>(), sFalse) && !args["path"]) {
        cerr << args.program << ": 'if --default FALSE' then  --path ('-p') option required" << endl;
        exit(1);
    }

    if (!in_array(args["default"].as<std::string>(), sFalse) && !in_array(args["default"].as<std::string>(), sTrue)) {
        cerr << args.program << ": '--default as FALSE or TRUE option required" << endl;
        exit(1);
    }
    if (!in_array(args["frequency"].as<std::string>(), sFalse) &&
        !in_array(args["frequency"].as<std::string>(), sTrue)) {
        cerr << args.program << ": '--frequency as FALSE or TRUE option required" << endl;
        exit(1);
    }

    if (in_array(args["default"].as<std::string>(), sFalse)) {
        cli.bDefaultNetwork = false;
        CheckPathExist(args["path"].as<std::string>());
        cli.sNetworkFile = args["path"].as<std::string>();
        cli.sNetworkName = args["network"].as<std::string>();
    }
    if (in_array(args["default"].as<std::string>(), sTrue)) {
        cli.bDefaultNetwork = true;
        cli.sNetworkFile = "";
        if (in_array(args["network"].as<std::string>(), sDefaultNetworks)) {
            cli.sNetworkName = args["network"].as<std::string>();
        } else {
            cerr << args.program << ": 'if --default TRUE' then --network can be only STRING10 or HumanNet1" << endl;
            exit(1);
        }
    }
    if (in_array(args["frequency"].as<std::string>(), sFalse)) {
        cli.bMutationFrequency = false;
    }
    if (in_array(args["frequency"].as<std::string>(), sTrue)) {
        cli.bMutationFrequency = true;
    }
    CheckPathExist(args["input"].as<std::string>());
    cli.sInputFile = args["input"].as<std::string>();
    cli.sOutputFile = args["output"].as<std::string>();

    return true;

}
