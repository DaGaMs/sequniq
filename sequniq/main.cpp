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

/* CHUNK is the size of the memory chunk used by the zlib routines. */

#define CHUNK 0x4000

/* The following macro calls a zlib routine and checks the return
 value. If the return value ("status") is not OK, it prints an error
 message and exits the program. Zlib's error statuses are all less
 than zero. */

#define CALL_ZLIB(x) {                                              \
    int status;                                                     \
    status = x;                                                     \
    if (status < 0) {                                               \
        fprintf (stderr,                                            \
                 "%s:%d: %s returned a bad status of %d.\n",        \
                 __FILE__, __LINE__, #x, status);                   \
        exit (EXIT_FAILURE);                                        \
    }                                                               \
}

/* These are parameters to deflateInit2. See
 http://zlib.net/manual.html for the exact meanings. */

#define windowBits 15
#define GZIP_ENCODING 16

static void strm_init (z_stream * strm)
{
    strm->zalloc = Z_NULL;
    strm->zfree  = Z_NULL;
    strm->opaque = Z_NULL;
    CALL_ZLIB (deflateInit2 (strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED,
                             windowBits | GZIP_ENCODING, 8,
                             Z_DEFAULT_STRATEGY));
}

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

void compress_to_stream (char *message, FILE *destination) {
    unsigned char out[CHUNK];
    z_stream strm;
    strm_init (& strm);
    strm.next_in = (unsigned char *) message;
    strm.avail_in = strlen (message);
    do {
        int have;
        strm.avail_out = CHUNK;
        strm.next_out = out;
        CALL_ZLIB (deflate (& strm, Z_FINISH));
        have = CHUNK - strm.avail_out;
        fwrite (out, sizeof (char), have, destination);
    }
    while (strm.avail_out == 0);
    deflateEnd (& strm);
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
    bool hasName = name != "";
    
    std::unordered_map<MurmurHash128, OffsetPair> hashtable;

    long lastOffset1 = 0;
    long lastOffset2 = 0;
    
    char *newChar = NULL;
    while ((l1 = kseq_read(seq1)) >= 0) {
        //calculate the hash function and print
        uint64_t *hashKey = new uint64_t[sizeof(uint64_t)*2];
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
        mh.size = 2;
        
        if (hashtable.count(mh)) {
            OffsetPair oldOP = hashtable.at(mh);
            if (oldOP.qual < op.qual)
                hashtable[mh] = op;
        } else {
            hashtable[mh] = op;
        }
    }
    free(newChar);
    
    FILE *output1 = stdout;
    FILE *output2 = NULL;
    
    if (hasName) {
        std::string outfileName;
        
        if (fp2)
            outfileName = name + "_1.fastq";
        else
            outfileName = outfileName + ".fastq";
        
        if (gzip)
            outfileName.append(".gz");

        output1 = fopen(outfileName.c_str(), "w");
        
        if (fp2)
            outfileName = name + "_2.fastq";
        
        if (gzip)
            outfileName.append(".gz");

        output2 = fopen(outfileName.c_str(), "w");
    }
    
    const int OUTBUFLEN = 262144;
    char *writeBuf1 = new char[sizeof(char)*OUTBUFLEN]; // 1MB output buffer
    char *writeBuf2 = NULL;
    if (fp2 && hasName) {
        writeBuf2 = new char[sizeof(char)*OUTBUFLEN]; // 1MB output buffer
    }
    char *strBuf =new char[sizeof(char) * 4096]; // 4000 chars line buf

    for(std::unordered_map<MurmurHash128, OffsetPair>::iterator iterator = hashtable.begin(); iterator != hashtable.end(); iterator++) {
        OffsetPair offsets = iterator->second;
        
        gzseek(fp1, offsets.offset1, SEEK_SET);
        seq1->last_char = 0;
        seq1->f->begin=0;
        seq1->f->end=0;
        seq1->f->is_eof=0;
        
        if ((l1 = kseq_read(seq1)) >= 0) {
            sprintf(strBuf, "@%s\n%s\n+\n%s\n", seq1->name.s, seq1->seq.s, seq1->qual.s);
        }
        
        if (strlen(writeBuf1) + strlen(strBuf) < OUTBUFLEN) {
            strcat(writeBuf1, strBuf);
        }
        else
        {
            if (gzip)
            {
                compress_to_stream(writeBuf1, output1);
            }
            else
            {
                fputs(writeBuf1, output1);
            }
            
            strcpy(writeBuf1, strBuf);
        }
        
        if (fp2)
        {
            gzseek(fp2, offsets.offset2, SEEK_SET);
            seq2->last_char = 0;
            seq2->f->begin=0;
            seq2->f->end=0;
            seq2->f->is_eof=0;
            if ((l2 = kseq_read(seq2)) >= 0) {
                sprintf(strBuf, "@%s\n%s\n+\n%s\n", seq2->name.s, seq2->seq.s, seq2->qual.s);
            }
            
            char *buf;
            if (hasName)
                buf = writeBuf2;
            else
                buf = writeBuf1;
            
            if (strlen(buf) + strlen(strBuf) < OUTBUFLEN) {
                strcat(buf, strBuf);
            }
            else
            {
                if (gzip)
                {
                    compress_to_stream(buf, hasName?output2:output1);
                }
                else
                {
                    fputs(buf, hasName?output2:output1);
                }
                
                strcpy(buf, strBuf);
            }
        }
    }
    
    if (strlen(writeBuf1))
    {
        if (gzip)
        {
            compress_to_stream(writeBuf1, output1);
        }
        else
        {
            fputs(writeBuf1, output1);
        }
    }
    
    if (fp2 && hasName)
    {
        if (strlen(writeBuf2))
        {
            if (gzip)
            {
                compress_to_stream(writeBuf2, output2);
            }
            else
            {
                fputs(writeBuf2, output2);
            }
        }
    }
    
    if (hasName) {
        fclose(output1);
        
        if (output2)
            fclose(output2);
    }
    
    free(writeBuf1);
    if (writeBuf2)
        free(writeBuf2);
    free(strBuf);
    kseq_destroy(seq1);
    kseq_destroy(seq2);
    gzclose(fp1);
    gzclose(fp2);
    return 0;
}


