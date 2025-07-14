# generate chunks for input fasta
# save gene name and corrdinates for each chunk in the header

import numpy as np
import os
from seekr.fasta_reader import Reader as seekrReader
from natsort import natsorted
import glob

# divide sequence into chunks around a target length without overlap
# length of seq divided by target length, the quotient is the number of chunks
# if remainder is greater than flex_length, add one more chunk
# if remainder is less than flex_length, the total num of chunk is the quotient
# divide the sequence equally into the total num of chunks, the remiander goes to the last chunk
def divide_seq(seq, target_length, flex_length): 
    seq_length = len(seq)
    quotient = seq_length // target_length
    remainder = seq_length % target_length
    if remainder > flex_length:
        num_chunk = quotient + 1
    else:
        num_chunk = quotient
    chunk_length = seq_length // num_chunk
    chunks = [seq[i:i+chunk_length] for i in range(0, seq_length, chunk_length)]
    # save the index range for each chunk
    # if last index is greater than seq_length, save seq_length instead
    chunk_index = [(i, np.min([i+chunk_length, seq_length])) for i in range(0, seq_length, chunk_length)]
    
    # append the last chunk to the second to the last chunk
    # if the last chunk length is less than chunk_length
    if len(chunks[-1]) < chunk_length:
        chunks[-2] = chunks[-2] + chunks[-1]
        # update chunk_index accordingly
        chunk_index[-2] = (chunk_index[-2][0], chunk_index[-1][1])
        del chunks[-1]
        del chunk_index[-1]
    return chunks, chunk_index


# read in fasta file and divide all seqs using divide_seq function
# write each divided chunks into a list

def divide_fasta_to_chunks(seqpath, target_length, flex_length): 
    seqs = seekrReader(seqpath).get_seqs()
    #headers=seekrReader(seqpath).get_headers()
    basename = os.path.basename(seqpath)
    name = os.path.splitext(basename)[0]
    # strip the "pos/neg_merged100" suffix for cleaner headers:
    prefix = name[:-14]
    div_seqs = []
    div_seq_headers = []
    for i, seq in enumerate(seqs):
        # only keep seqs longer than flex_length
        #if len(seq) > flex_length:
        # keep all seqs
        if len(seq) > 30:
            # if seq is longer than target_length, divide it into chunks
            if len(seq) > target_length:
                chunks, chunk_index = divide_seq(seq, target_length, flex_length)
                div_seqs.append(chunks)
                # save the header for all chunks
                # chunk_headers=[headers[i] +'_' + str(index[0]) + '_' + str(index[1]) for index in chunk_index]
                chunk_headers = [
                    '>' + prefix +'_p' + str(i+1) + 'chk' + str(j + 1)
                    for j in range(len(chunk_index))
                ]
                div_seq_headers.append(chunk_headers)
            # if seq is shorter than target_length, append it to the list
            else:
                div_seqs.append([seq])
                # add seq length to the header
                div_seq_headers.append(['>' + prefix + '_p'+ str(i+1) + 'chk1'])
    # flatten the list of lists into a list
    div_seqs = [item for sublist in div_seqs for item in sublist]
    div_seq_headers = [item for sublist in div_seq_headers for item in sublist]
    return div_seqs, div_seq_headers


# save all element of a list of seqs into one fasta file
# close the files after written to avoid too many open files error

def save_seqs_to_fasta(seqs, seq_names, seqpath):
    seqfile = open(seqpath, 'w')
    for i, seq in enumerate(seqs):
        # save fasta name from corresponding element in seq_names
        seqfile.write(seq_names[i] + '\n' + seq + '\n')
    seqfile.close()


fa_files = natsorted(glob.glob("*_merged100.fa"))

for filepath in fa_files:
    basename = os.path.basename(filepath)
    name = os.path.splitext(basename)[0]

    div_seqs, div_seq_headers = divide_fasta_to_chunks(filepath, 200, 150)
    # save  into a fasta file
    save_seqs_to_fasta(div_seqs, div_seq_headers, f"{name}_200cut.fa")



# change how the chunks are named
def divide_fasta_to_chunks(seqpath, target_length, flex_length): 
    seqs = seekrReader(seqpath).get_seqs()
    #headers=seekrReader(seqpath).get_headers()
    basename = os.path.basename(seqpath)
    name = os.path.splitext(basename)[0]
    # only keep the first 4 letter:
    prefix = name[:4]+'_null'
    div_seqs = []
    div_seq_headers = []
    for i, seq in enumerate(seqs):
        # only keep seqs longer than flex_length
        #if len(seq) > flex_length:
        # keep all seqs
        if len(seq) > 30:
            # if seq is longer than target_length, divide it into chunks
            if len(seq) > target_length:
                chunks, chunk_index = divide_seq(seq, target_length, flex_length)
                div_seqs.append(chunks)
                # save the header for all chunks
                # chunk_headers=[headers[i] +'_' + str(index[0]) + '_' + str(index[1]) for index in chunk_index]
                chunk_headers = [
                    '>' + prefix +'_p' + str(i+1) + 'chk' + str(j + 1)
                    for j in range(len(chunk_index))
                ]
                div_seq_headers.append(chunk_headers)
            # if seq is shorter than target_length, append it to the list
            else:
                div_seqs.append([seq])
                # add seq length to the header
                div_seq_headers.append(['>' + prefix + '_p'+ str(i+1) + 'chk1'])
    # flatten the list of lists into a list
    div_seqs = [item for sublist in div_seqs for item in sublist]
    div_seq_headers = [item for sublist in div_seq_headers for item in sublist]
    return div_seqs, div_seq_headers


    
div_seqs, div_seq_headers = divide_fasta_to_chunks('Xist_null.fa', 200, 150)
save_seqs_to_fasta(div_seqs, div_seq_headers, f"Xist_null_200cut.fa")


div_seqs, div_seq_headers = divide_fasta_to_chunks('Airn_null.fa', 200, 150)
save_seqs_to_fasta(div_seqs, div_seq_headers, f"Airn_null_200cut.fa")


div_seqs, div_seq_headers = divide_fasta_to_chunks('Kot1_null.fa', 200, 150)
save_seqs_to_fasta(div_seqs, div_seq_headers, f"Kot1_null_200cut.fa")

