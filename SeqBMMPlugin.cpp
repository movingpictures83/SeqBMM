#include "PluginManager.h"
#include <stdio.h>
#include <stdlib.h>
#include "SeqBMMPlugin.h"

void SeqBMMPlugin::input(std::string file) {
 inputfile = file;
 std::ifstream ifile(inputfile.c_str(), std::ios::in);
 while (!ifile.eof()) {
   std::string key, value;
   ifile >> key;
   ifile >> value;
   parameters[key] = value;
 }
}

void SeqBMMPlugin::run() {

}


void SeqBMMPlugin::output(std::string file) {
   std::string command = "export OLDPATH=${PYTHONPATH}; ";
   command += "export PYTHONPATH=/usr/local/lib64/python3.9/site-packages/:${PYTHONPATH}; ";
   command += "python3.9 plugins/SeqBMM/runSeqBMM.py ";
   command += PluginManager::addPrefix(parameters["sqldatabase"]) + " ";
   command += parameters["pdbinput"] + " ";
   command += parameters["span"] + " ";
   command += PluginManager::addPrefix(parameters["pdbdatabase"]) + " ";
   command += PluginManager::addPrefix(parameters["rasafile"]) + " ";
   command += PluginManager::addPrefix(parameters["seqnumsfile"]) + " ";
   command += PluginManager::addPrefix(parameters["seqsolvfile"]) + " ";
   command += PluginManager::addPrefix(parameters["csvfile"]) + " ";
   command += file + "; ";
   command += "export PYTHONPATH=${OLDPATH}";
 std::cout << command << std::endl;

 system(command.c_str());
}

PluginProxy<SeqBMMPlugin> SeqBMMPluginProxy = PluginProxy<SeqBMMPlugin>("SeqBMM", PluginManager::getInstance());
