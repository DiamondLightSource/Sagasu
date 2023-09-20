. /etc/profile.d/modules.sh

module load phenix

phenix.emma --symmetry="$1" "$2" "$3"
