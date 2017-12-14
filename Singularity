Bootstrap: docker
From: debian:jessie

%post
#### install system dependencies
apt-get update
apt-get clean
apt-get install --no-install-recommends -qy \
                  cmake \
                  g++-4.9 \
                  git \
                  libboost-all-dev \
                  make \
                  python3 \
                  wget

#### install bpp
cd /opt/

echo "deb http://biopp.univ-montp2.fr/repos/apt/ Trusty main" >> /etc/apt/sources.list;
wget http://biopp.univ-montp2.fr/repos/apt/conf/biopp.gpg.key
apt-key add biopp.gpg.key
apt-get update
apt-get install -qy libbpp-phyl-dev

#### compile and install ALE
git clone git://github.com/ssolo/ALE /usr/local/ALE
mkdir /usr/local/ALE/build

cd /usr/local/ALE/build

echo "export LD_LIBRARY_PATH=/usr/local/lib/" >> $SINGULARITY_ENVIRONMENT

cmake ..  -DCMAKE_CXX_COMPILER=/usr/bin/g++-4.9 && make -j 4

for binary in /usr/local/ALE/build/bin/*; do ln -s $binary /usr/local/bin/; done

%labels
Maintainer	max-emil.schon@icm.uu.se
