#ifndef SEQBMMPLUGIN_H
#define SEQBMMPLUGIN_H

#include "Plugin.h"
#include "PluginProxy.h"
#include <string>
#include <vector>

class SeqBMMPlugin : public Plugin
{
public: 
 std::string toString() {return "SeqBMM";}
 void input(std::string file);
 void run();
 void output(std::string file);

private: 
 std::string inputfile;
 std::string outputfile;
 std::map<std::string, std::string> parameters;

};

#endif
