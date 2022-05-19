#install doxygen using the command bellow
sudo apt-get install doxygen

# create doxygen config file using the command bellow
doxygen -s -g myConfigFileName

Edit your config file using your favorite editor. 
The important fields to be updated are:

PROJECT_NAME           = "MyProject"

OUTPUT_DIRECTORY       =../workspace/myFavoriteDirectory
INPUT                  =../workspace/MyFavProjectLink/
RECURSIVE              = YES
GENERATE_HTML          = YES

#The config file is written in very human readable format, so make any other changes you would like.

# To generate doxygen documentation run the bellow command
doxygen myConfigFileName

