#include <zlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <algorithm>
#include <unordered_map>

#include "kseq.h"
#include "MurmurHash3.h"
#include <tclap/CmdLine.h>

#if __x86_64__
/* 64-bit */
#define murmur(k, l, s, o) MurmurHash3_x64_128(k, l, s, o)
#else
#define murmur(k, l, s, o) MurmurHash3_x86_128(k, l, s, o)
#endif

//using namespace std;

KSEQ_INIT(gzFile, gzread)

std::hash<uint64_t> longHash;

struct MurmurHash128
{
    uint64_t* p;
    size_t size;
    
    bool operator==(const MurmurHash128 &other) const
    {
        return size == other.size && memcmp(p, other.p, size) == 0;
    }
};

namespace std {
    template <>
    struct hash<MurmurHash128>
    {
        typedef MurmurHash128 argument_type;
        typedef std::size_t result_type;
        result_type operator()(argument_type const& dr) const
        {
            result_type h = 0;
            for (uint64_t* p = dr.p; p != dr.p + dr.size; ++p)
                h ^= longHash(*p);
            return h;
        }
    };
}

struct OffsetPair
{
    long offset1;
    long offset2;
    int qual;
};

int calculate_score (const char *scoreString) {
    int result = 0;
    while (*scoreString){
        result += *scoreString - 33;
        scoreString++;
    }
    return result;
}

int main(int argc, char *argv[])
{
    std::string name;
    bool gzip;

    std::string input1file, input2file;
    
    try {
        TCLAP::CmdLine cmd("sequniq removes identical reads from FastQ files, retaining only the highest score string.", ' ', "0.1");
        TCLAP::ValueArg<std::string> nameArg("p","prefix","File name prefix for output",false,"","prefix");
        cmd.add( nameArg );
        
        TCLAP::SwitchArg gzipSwitch("z","gzip","Compress output", false);
        cmd.add( gzipSwitch );
        
        TCLAP::UnlabeledValueArg<std::string> input1arg("file1.fq[.gz]", "FastQ file (optionally gzip compressed) to be filtered", true, "", "file1.fq[.gz]", cmd);
        TCLAP::UnlabeledValueArg<std::string> input2arg("file2.fq[.gz]", "FastQ file (optionally gzip compressed) with paired reads to file 1", false, "", "file2.fq[.gz]", cmd);

        // Parse the argv array.
        cmd.parse( argc, argv );
        
        // Get the value parsed by each arg.
        name = nameArg.getValue();
        gzip = gzipSwitch.getValue();
        input1file = input1arg.getValue();
        input2file = input2arg.getValue();
    } catch (TCLAP::ArgException &e)  // catch any exceptions
    { std::cerr << "ERROR: " << e.error() << " for arg " << e.argId() << std::endl; }

    uint32_t seed = rand(); // random hash seed
    
    gzFile fp1, fp2;
    kseq_t *seq1, *seq2;
    
    int l1, l2;

    fp1 = gzopen(input1file.c_str(), "r");
    seq1 = kseq_init(fp1);
    
    if (input2file != "") {
        fp2 = gzopen(input2file.c_str(), "r");
        seq2 = kseq_init(fp2);
    }
    else
    {
        fp2 = NULL;
        seq2 = NULL;
    }
    
    std::unordered_map<MurmurHash128, OffsetPair> hashtable;

    long lastOffset1 = 0;
    long lastOffset2 = 0;
    
    char *newChar = NULL;
    while ((l1 = kseq_read(seq1)) >= 0) {
        //calculate the hash function and print
        uint64_t *hashKey = new uint64_t[2];
        MurmurHash128 mh;
        OffsetPair op;
        
        if (fp2)
        {
            if ((l2 = kseq_read(seq2)) >= 0)
            {
                if (sizeof(newChar) < (strlen(seq1->seq.s) + strlen(seq2->seq.s) + 1))
                {
                    free(newChar);
                    newChar = NULL;
                }

                if (!newChar) {
                    newChar = new char[strlen(seq1->seq.s) + strlen(seq2->seq.s) + 1];
                }

                strcpy(newChar, seq1->seq.s);
                strcat(newChar, seq2->seq.s);
                murmur(newChar, (int)strlen(newChar), seed, hashKey);
            }
            else
            {
                fprintf(stderr, "ERROR: paired-end files have different length");
                return 2;
            }
            
            op.offset1 = lastOffset1;
            op.offset2 = lastOffset2;

            lastOffset1 = (gztell(fp1) - seq1->f->end) + seq1->f->begin;
            lastOffset2 = (gztell(fp2) - seq2->f->end) + seq2->f->begin;
            
            op.qual = op.qual = calculate_score(seq1->qual.s) + calculate_score(seq2->qual.s);
        }
        else
        {
            murmur(seq1->seq.s, (int)seq1->seq.l, seed, hashKey);
            op.offset1 = lastOffset1;
            lastOffset1 = (gztell(fp1) - seq1->f->end) + seq1->f->begin;
            
            op.qual = calculate_score(seq1->qual.s);
        }
        
        mh.p = hashKey;
        mh.size = 128;
        
        if (hashtable.count(mh)) {
            OffsetPair oldOP = hashtable.at(mh);
            if (oldOP.qual < op.qual)
                hashtable[mh] = op;
        } else {
            hashtable[mh] = op;
        }
    }
    free(newChar);
    
    for(std::unordered_map<MurmurHash128, OffsetPair>::iterator iterator = hashtable.begin(); iterator != hashtable.end(); iterator++) {
        OffsetPair offsets = iterator->second;
        
        gzseek(fp1, offsets.offset1, SEEK_SET);
        seq1->last_char = 0;
        seq1->f->begin=0;
        seq1->f->end=0;
        seq1->f->is_eof=0;
        if ((l1 = kseq_read(seq1)) >= 0) {
            printf("@%s\n%s\n+\n%s\n", seq1->name.s, seq1->seq.s, seq1->qual.s);
        }
        
        if (fp2)
        {
            gzseek(fp2, offsets.offset2, SEEK_SET);
            seq2->last_char = 0;
            seq2->f->begin=0;
            seq2->f->end=0;
            seq2->f->is_eof=0;
            if ((l2 = kseq_read(seq2)) >= 0) {
                printf("@%s\n%s\n+\n%s\n", seq2->name.s, seq2->seq.s, seq2->qual.s);
            }
        }
    }
    
    kseq_destroy(seq1);
    kseq_destroy(seq2);
    gzclose(fp1);
    gzclose(fp2);
    return 0;
}
