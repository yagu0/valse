#!/bin/sh

#initialize submodules, set-up .git/config and .gitattributes, and pre-push hook
git submodule init && git submodule update --merge
#filter for Jupyter
python .nbstripout/nbstripout.py --install --attributes .gitattributes
#filter for git-fat [TODO: idempotent...]
printf '*.pdf filter=fat\n*.tar.xz filter=fat\n*.png filter=fat\n*.jpg filter=fat\n*.ps filter=fat\n'  >> .gitattributes
#pre-push hook: git fat push, submodules update
printf '#!/bin/sh\n./.git-fat/git-fat pull\n./.git-fat/git-fat push\ngit submodule update --merge\n' > .git/hooks/pre-push
chmod 755 .git/hooks/pre-push
#.gitfat file with remote on gitfat@auder.net
printf '[rsync]\nremote = gitfat@auder.net:~/files/valse\n' > .gitfat
#manual git-fat init: with relative path to binary
#1] remove filter if exists http://stackoverflow.com/questions/12179437/replace-3-lines-with-another-line-sed-syntax
sed -i '1N;$!N;s/\[filter "fat"\]\n.*\n.*//;P;D' .git/config
#2] place new filter
printf '[filter "fat"]\n\tclean = ./.git-fat/git-fat filter-clean\n\tsmudge = ./.git-fat/git-fat filter-smudge\n' >> .git/config
