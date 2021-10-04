#!/usr/bin/env bash

if ! [ -f araTha1/araTha1.fa.sizes ]; then
	mkdir -p araTha1
	cd araTha1
	wget https://genome-test.gi.ucsc.edu/~hiram/hubs/Plants/araTha1/araTha1.2bit
	twoBitToFa araTha1.2bit araTha1.fa
	faidx -x araTha1.fa
	faidx araTha1.fa -i chromsizes > araTha1.fa.sizes
	rm araTha1.2bit chr*.fa
	cd ..
fi

if ! [ -f ce10/ce10.fa.sizes ]; then
	genomepy install -p UCSC -g ./ -r "^chr[XVIM]{1,3}$" -t 8 -f ce10
fi

if ! [ -f ce11/ce11.fa.sizes ]; then
	genomepy install -p UCSC -g ./ -r "^chr[XVIM]{1,3}$" -t 8 -f ce11
fi

if ! [ -f ci3/ci3.fa.sizes ]; then
	genomepy install -p UCSC -g ./ -r "^chr[\d]{1,2}$" -t 8 -f ci3
fi

if ! [ -f danRer11/danRer11.fa.sizes ]; then
	genomepy install -p UCSC -g ./ -r "^chr[\dXYM]{1,2}$" -t 8 -f danRer11
fi

if ! [ -f dm6/dm6.fa.sizes ]; then
	genomepy install -p UCSC -g ./ -r "^chr[234XYLR]{1,2}$" -t 8 -f dm6
fi

if ! [ -f hg19/hg19.fa.sizes ]; then
	genomepy install -p UCSC -g ./ -r "^chr[\dXYM]{1,2}$" -t 8 -f hg19
fi

if ! [ -f hg38/hg38.fa.sizes ]; then
	genomepy install -p UCSC -g ./ -r "^chr[\dXYM]{1,2}$" -t 8 -f hg38
fi

if ! [ -f mm10/mm10.fa.sizes ]; then
	genomepy install -p UCSC -g ./ -r "^chr[\dXYM]{1,2}$" -t 8 -f mm10
fi

if ! [ -f mm39/mm39.fa.sizes ]; then
	genomepy install -p UCSC -g ./ -r "^chr[\dXYM]{1,2}$" -t 8 -f mm39
fi

if ! [ -f sacCer3/sacCer3.fa.sizes ]; then
	genomepy install -p UCSC -g ./ -r "^chr[XVIM]{1,4}$" -t 8 -f sacCer3
fi
