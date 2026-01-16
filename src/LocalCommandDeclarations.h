#ifndef LOCALCOMMANDDECLARATIONS_H
#define LOCALCOMMANDDECLARATIONS_H

#include "Command.h"

//extern int easypredict(int argc, const char **argv, const Command& command);
extern int createsetdb(int argc, const char **argv, const Command& command);
extern int gff2db(int argc, const char **argv, const Command& command);
extern int aa2foldseek(int argc, const char **argv, const Command& command);
extern int clustersearch(int argc, const char **argv, const Command& command);
extern int clusterhits(int argc, const char **argv, const Command& command);
extern int besthitbyset(int argc, const char **argv, const Command &command);
extern int combinehits(int argc, const char **argv, const Command &command);
extern int summarizeresults(int argc, const char **argv, const Command &command);
extern int clusterdb(int argc, const char **argv, const Command &command);
#endif
