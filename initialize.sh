#!/bin/sh

#initialize submodules, set-up .git/config and .gitattributes, and pre-push hook
git submodule init && git submodule update --merge

#filter for git-fat
printf \
'*.pdf filter=fat
*.tar.xz filter=fat
*.png filter=fat
*.jpg filter=fat
*.ps filter=fat\n' > .gitattributes

#filter for Jupyter
python .nbstripout/nbstripout.py --install --attributes .gitattributes

#pre-commit and pre-push hooks: indentation, git fat push, submodules update
cp hooks/* .git/hooks/

#install formatR
echo 'if (! "formatR" %in% rownames(installed.packages()))
	install.packages("formatR",repos="https://cloud.r-project.org")' | R --slave

#.gitfat file with remote on gitfat@auder.net
printf '[rsync]\nremote = gitfat@auder.net:~/files/valse\n' > .gitfat

#manual git-fat init: with relative path to binary
#1] remove filter if exists http://stackoverflow.com/a/12179641/4640434
sed -i '1N;$!N;s/\[filter "fat"\]\n.*\n.*//;P;D' .git/config
#2] place new filter
printf \
'[filter "fat"]
	clean = ./.git-fat/git-fat filter-clean
	smudge = ./.git-fat/git-fat filter-smudge\n' >> .git/config
