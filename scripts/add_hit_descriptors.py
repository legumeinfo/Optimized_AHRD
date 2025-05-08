#!/usr/bin/env python3
import sys
import csv
import re
import os

def main():
    
    ahrd_output_path = sys.argv[1]
    go_lookup_path = sys.argv[2]
    out_dir = sys.argv[3]
    output_path = os.path.join(out_dir, "ahrd_with_all_descriptions.tsv")
    
    go_descriptions = {}
    with open(go_lookup_path, 'r', newline='') as lookup_file:
        reader = csv.reader(lookup_file, delimiter='\t')
        # skips header automaticaly
        next(reader)
        for row in reader:
            # lookup file should have 2 cols
            if len(row) >= 2:
                go_id = row[0]
                description = row[1]
                go_descriptions[go_id] = description
    
    # pulls GO terms via regex
    go_term_pattern = re.compile(r'GO:\d{7}')
    
    with open(ahrd_output_path, 'r', newline='') as input_file, \
         open(output_path, 'w', newline='') as output_file:
        
        reader = csv.reader(input_file, delimiter='\t')
        writer = csv.writer(output_file, delimiter='\t')
        
        # skips header
        header = next(reader)
        writer.writerow(header)
        # Process data rows
        for row in reader:
            
            if len(row) >= 6:  # makes sure at least 6 cols
                go_terms_field = row[5] if len(row) > 5 else ""
                
                if go_terms_field:
                    # splits by comma followed by optional space
                    terms = [term.strip() for term in go_terms_field.split(',')]
                    enhanced_terms = []
                    
                    for term in terms:
                        # checks it's a valid GO term
                        match = go_term_pattern.search(term)
                        if match:
                            go_id = match.group(0)
                            if go_id in go_descriptions:
                                enhanced_term = f"{term} ({go_descriptions[go_id]})"
                                enhanced_terms.append(enhanced_term)
                            #else:
                                #enhanced_terms.append(term)
                        else:
                            enhanced_terms.append(term)
                    
                    # modifies the column
                    row[5] = ', '.join(enhanced_terms)
            
            writer.writerow(row)

if __name__ == "__main__":
    main()
