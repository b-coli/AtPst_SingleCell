sudo apt-get update
sudo apt-get install -y \
  build-essential \
  libseccomp-dev \
  libglib2.0-dev \
  pkg-config \
  squashfs-tools \
  cryptsetup runc \
  pigz

export VERSION=1.18.4 OS=linux ARCH=amd64

wget -O /tmp/go${VERSION}.${OS}-${ARCH}.tar.gz \
  https://dl.google.com/go/go${VERSION}.${OS}-${ARCH}.tar.gz

sudo tar -C /usr/local -xzf /tmp/go${VERSION}.${OS}-${ARCH}.tar.gz

echo 'export PATH=$PATH:/usr/local/go/bin' >> ~/.bashrc
source ~/.bashrc

git clone --recurse-submodules https://github.com/sylabs/singularity.git

cd singularity

git checkout --recurse-submodules v3.10.2

./mconfig
make -C builddir
sudo make -C builddir install

cd ~/

singularity pull docker://bcoli/atpst_singlecell_r:1.0
singularity pull docker://dceoy/gffread:latest
singularity pull docker://dceoy/trim_galore:latest
singularity pull docker://staphb/fastqc:0.11.9
singularity pull docker://cumulusprod/cellranger:6.0.1
singularity pull docker://mparikhbroad/velocyto:1.0.1
singularity pull docker://ncbi/sra-tools:3.0.0

