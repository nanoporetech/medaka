#ifndef _MEDAKA_FASTRLE_H
#define _MEDAKA_FASTRLE_H

/**  Read a fasta/q and writes quality-RLE fastq file.
 *
 *  @param fname name of file to be rle-compressed.
 *  @param block_size maximum RLE compressed size, homopolymers larger than this value will be split in blocks.
 *  @returns 
 *
 */
size_t fastrle(char *fname, size_t block_size);
#endif
