#! /usr/bin/python3
#Authors: Patrick Roelli, Kirchner Benedict
from itertools import product
from collections import defaultdict
import argparse
from argparse import RawTextHelpFormatter
from anytree import Node, RenderTree, PreOrderIter, LevelOrderIter
from anytree.importer import JsonImporter
import copy
import progressbar
import sys

def get_args():
    """Get args."""
    parser = argparse.ArgumentParser(prog='CITE Seq Count',
        description="""This creates a reference new isomir reference.""", formatter_class=RawTextHelpFormatter)
    inputs = parser.add_argument_group('Inputs')
    inputs.add_argument('-H', '--hairpin',
                        action='store',
                        help='The path to hairpin file',
                        dest='hairpin_path',
                        required=True)
    inputs.add_argument('-M', '--mature',
                        help='The path mature ref',
                        dest='mature_path',
                        required=True)
    inputs.add_argument('-3T', '--3_prime_template',
                        help='The 3prime end tree template as json',
                        dest='three_prime_template_file',
                        required=True)
    inputs.add_argument('-5T', '--5_prime_template',
                        help='The 5 prime end tree template as json',
                        dest='five_prime_template_file',
                        required=True)
    inputs.add_argument('-OUT', '--output_file',
                        help='File path for isomir reference',
                        dest='out_path',
                        required=True)
    inputs.add_argument('-MOUT', '--multireads',
                        help='File path for multireads',
                        dest='multireads_path',
                        required=True)
    return parser


def match_one_end(sequence, reference_dict, end, hairpin_reference,max_add=3):
    """Match one end to hairpin reference"""
    for node in LevelOrderIter(reference_dict['{}_prime_tree'.format(end)], maxlevel=1):
        node.sequence = sequence
    for node in PreOrderIter(reference_dict['{}_prime_tree'.format(end)], maxlevel=max_add+1):
        if(not node.is_root):
            if(end==5):
                node.sequence = ''.join([node.base,node.parent.sequence])
            elif(end==3):
                node.sequence = ''.join([node.parent.sequence,node.base])
            for hairpin in hairpin_reference:
                if (node.sequence in hairpin):
                    node.is_template = True
                else:
                    #TODO: Find how to stop iterating if false and jump to next branch
                    continue
    return(reference_dict)


def generate_half_string(sequence, reference, part, type):
    """Generate half of the name to add"""
    string_list = []
    if(type == 'Add'):
        for node in PreOrderIter(reference['{}_prime_tree'.format(part)]):
            if (not node.is_root):
                string_list.append('{}_{}_{}_{}'.format(
                    '{}{}'.format('Add',part),
                    len(node.name),
                    node.name,
                    ['nontemplate','template'][node.is_template]))
    elif(type == 'Can'):
        string_list.append('{}_{}_{}_{}'.format(
            '{}{}'.format('Can',part),
            0,
            'X',
            'template'))
    elif(type == 'Trim'):
        for number in range(1,7,1):
            string_list.append('{}_-{}_{}_{}'.format(
                    '{}{}'.format('Trim',part),
                    number,
                    [sequence[-number:],sequence[:number]][part==5],
                    'template'))
    return(string_list)

    
def modify_sequence(sequence, mod5, mod3):
    """Modify the sequence based on the name"""
    new_seq = sequence
    mod5 = mod5.split('_')
    mod3 = mod3.split('_')
    
    if(int(mod5[1])>0):
        new_seq = mod5[2] + new_seq
    if(int(mod5[1])<0):
        new_seq = new_seq[-int(mod5[1]):]
    
    if(int(mod3[1])>0):
        new_seq = new_seq + mod3[2]
    if(int(mod3[1])<0):
        new_seq = new_seq[:int(mod3[1])]
    return(new_seq)

def generate_string_general(sequence, reference, new_reference):
    """Create all the combinations of names"""
    types = ['Add','Trim','Can']
    five_prime_parts=[]
    three_prime_parts=[]
    for type in types:
        five_prime_parts.extend(generate_half_string(sequence=sequence, reference=reference[sequence], part=5,type=type))
        three_prime_parts.extend(generate_half_string(sequence=sequence, reference=reference[sequence], part=3,type=type))

    all_comb = product(five_prime_parts, three_prime_parts)

    for comb in all_comb:
        seq = modify_sequence(sequence, comb[0],comb[1])
        new_name = '{}_{}_{}'.format(reference[sequence]['name'],comb[0],comb[1])
        if seq in new_reference:
            new_reference[seq].add(new_name)
        else:
            new_reference[seq] = set()
            new_reference[seq].add(new_name)
    return(new_reference)


def order_multireads(names):
    """Order the multireads"""
    names_list = list(names)
    names_split=[]
    for name in names_list:
        name_split = name.split('_')
        names_split.append(name_split)
    #This sorts a nested list. Uses reverse=True for alphabetical order
    #Uses -abs(int(X)) because we want it ordered normally, hence the minus sign
    names_split = sorted(names_split, key=lambda x: (-abs(int(x[2])), -abs(int(x[6])), x[1], x[5]), reverse=True)
    ordered_names = []
    for name in names_split:
        ordered_names.append('_'.join(name))
    return(ordered_names)


def main():
    parser = get_args()
    if not sys.argv[1:]:
        sys.exit(parser.print_help())
    #Load args
    args = parser.parse_args()
    n = 0
    importer = JsonImporter()
    with open(args.three_prime_template_file) as json_data:
        three_prime_template = importer.read(json_data)
    with open(args.five_prime_template_file) as json_data:
        five_prime_template = importer.read(json_data)
    hairpins_list=set()
    reference = defaultdict()
    with open(args.hairpin_path, 'rt') as hairpins, open(args.mature_path, 'rt') as mature: 
        #Load hairpins ref
        all_lines = hairpins.read()
        hairpins = all_lines.split('>')
        del hairpins[0]
        for line in hairpins:
            line = line.split('\n',1)
            hairpins_list.add(line[1].replace('\n',''))
        #Load mature RNA
        for line in mature:
            line2 = next(mature).strip()
            name = line.strip()
            if line2 in reference:
                reference[line2]['name'] = '{}={}'.format(reference[line2]['name'],name.lstrip('>')) 
            else:
                reference[line2]=defaultdict()
                reference[line2]['name'] = name
            reference[line2]['5_prime_tree'] = copy.deepcopy(five_prime_template)
            reference[line2]['3_prime_tree'] = copy.deepcopy(three_prime_template)
        num_seqs = len(reference.keys())
        widgets = [progressbar.Percentage(), progressbar.Bar()]
        print('Mapping isomir modifications to reference.')
        count=0    
        bar = progressbar.ProgressBar(widgets=widgets, max_value=num_seqs).start()
        for sequence in reference.keys():
            match_one_end(sequence=sequence,reference_dict=reference[sequence],end=5,hairpin_reference=hairpins_list)
            match_one_end(sequence=sequence,reference_dict=reference[sequence],end=3,hairpin_reference=hairpins_list)
            count+=1
            bar.update(count)
        print()
        print('Generating new reference.')
        count = 0
        new_reference = defaultdict()
        bar = progressbar.ProgressBar(widgets=widgets, max_value=num_seqs).start()
        for sequence in reference.keys():
            new_reference = generate_string_general(sequence, reference, new_reference)
            count+=1
            bar.update(count)
        print()
        print('Checking for multireads.')
        with open(args.out_path,'w') as out_file:
            with open(args.multireads_path,'w') as multireads_out:
                for sequence in new_reference:
                    if(len(new_reference[sequence])>1):
                        ordered_names = order_multireads(names=new_reference[sequence])
                        multireads_out.write('{}\n{}\n\n'.format('\n'.join(ordered_names), sequence))
                        out_file.write('{}_multi\n{}\n'.format(ordered_names[0],sequence))
                    else:
                        out_file.write('{}_unique\n{}\n'.format(next(iter(new_reference[sequence])),sequence))
        print('Isomir reference created.')

if __name__ == '__main__':
    main()