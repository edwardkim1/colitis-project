#!/bin/bash
#sshfs_setup.sh

echo "Unmount disk?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) echo "Disk to unmount:"; read pathToDir;diskutil unmount force $pathToDir;hdiutil eject -force $pathToDir; break;;
        No ) break;;
    esac
done

echo "Mount 1.projects?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) echo "Directory to mount:"; read toMount;sshfs -p 22 esk17@rcapps5.dfci.harvard.edu:/mnt/beegfs/home/esk17/1.projects/ $toMount; break;;
        No ) exit;;
    esac
done

