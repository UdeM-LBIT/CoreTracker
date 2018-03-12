# $1 = new version
# $2 = date
# update new version
tag="$1"
if [ -z "$tag" ]
then
	tag=`grep -Po  "\_\_version\_\_\s=\s\K.*" coretracker/__init__.py | tr -d "'" |  awk -F. -v OFS=. 'NF==1{print ++$NF}; NF>1{if(length($NF+1)>length($NF))$(NF-1)++; $NF=sprintf("%0*d", length($NF), ($NF+1)%(10^length($NF))); print}'`
fi

echo "New version $tag"
sed -i -e "s/\_\_version\_\_.*/\_\_version\_\_ = '$tag'/g" coretracker/__init__.py

if ! [ -z "$2" ]
then
	sed -i -e "s/date =.*/date = '$2'/g" coretracker/__init__.py
fi
# update conda build version
sed -i -e "s/{\% set version =.*\%}/{% set version = \"$tag\" %}/g" conda_build/meta.yaml

git add coretracker/__init__.py conda_build/meta.yaml
git commit -m "bumb version"
git push origin master
# now create a new git tag 
git tag $tag
git push origin $tag
python setup.py sdist
twine upload "dist/CoreTracker-""$tag"".tar.gz"
cd conda_build && conda build .
cd ..