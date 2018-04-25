#!/bin/bash
#creates json list for trees to chose from

#begin json list with left bracket
echo -n "[" > contents.json          
#look at current dir to find all newick files
#change ls *rooting.0 to directory name if collecting all files from a folder disregarding extension -- current usage looks for files within cwd
for newick in $(ls Phylotree_files); do
    echo -e -n "\"$newick\"," >> contents.json #append name of file along with a comma to json list
done
#close json list with right bracket
echo -n "]" >> contents.json
#replace instance of ,] with ] -> aka the last file added
sed -i -e 's/,]/]/g' contents.json

#To view viz uncomment the two lines below
#open http://127.0.0.1:8000/
#python -m SimpleHTTPServer
