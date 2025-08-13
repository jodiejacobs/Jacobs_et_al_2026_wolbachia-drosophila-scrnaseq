#!/usr/bin/env python3

# Genus colors
genera = [
    ('Tatumella', '#89a731'),
    ('Serratia', '#ad9c31'),
    ('Huaxiibacter', '#7e97f4'),
    ('Gibbsiella', '#34ae92'),
    ('Shewanella', '#35aca4'),
    ('Silvania', '#f56ab4'),
    ('Edwardsiella', '#ec7e32'),
    ('Citrobacter', '#f66f94'),
    ('Pectobacterium', '#f7754f'),
    ('Vibrio', '#e06cf4'),
    ('Enterobacter', '#c29431'),
    ('Gilliamella', '#5caf31'),
    ('Wolbachia', '#36abaf'),
    ('Bartonella', '#ae88f4'),
    ('Ehrlichia', '#c27ff4'),
    ('Anaplasma', '#34ae8d'),
    ('Neorickettsia', '#38a8c5'),
    ('Sphingomonas', '#90a531'),
    ('Acetobacter', '#80a831'),
    ('Brachybacterium', '#d48c31'),
    ('Thermoanaerobacterium', '#39a6d7'),
    ('Commensalibacter', '#96a331'),
    ('Novosphingobium', '#9890f4'),
    ('Orbus', '#b79931'),
    ('Thermohydrogenium', '#37aab7'),
    ('Pseudomonas', '#f563d4'),
    ('Gemella', '#37aabb'),
    ('Corynebacterium', '#d673f4'),
    ('Porticoccus', '#f77462'),
    ('Lachnoanaerobaculum', '#32b071'),
    ('Utexia', '#36abab'),
    ('Luteimonas', '#35ada0'),
    ('Cutibacterium', '#34ad97'),
    ('Neisseria', '#6e9af4'),
    ('Methylobacterium', '#eb63f4'),
    ('Thiomonas', '#8c93f4'),
    ('Ancylobacter', '#38a7d0'),
    ('Streptomyces', '#c89231'),
    ('Rhodococcus', '#3ca1f4'),
    ('Pedococcus', '#f567c4'),
    ('Lapillicoccus', '#31b256'),
    ('Dermacoccus', '#f66d9c'),
    ('Hyphomicrobium', '#f461dd'),
    ('Parasphingopyxis', '#47b131'),
    ('Pseudotabrizicola', '#35aca7'),
    ('Cypionkella', '#e38431'),
    ('Microbacterium', '#f565cc'),
    ('Acinetobacter', '#32b165'),
    ('Stenotrophomonas', '#a89e31'),
    ('Pseudoxanthomonas', '#3ba3e8'),
    ('Xanthomonas', '#36abb3'),
    ('Lysobacter', '#f66ca5'),
    ('Nitrosomonas', '#cc79f4'),
    ('Moraxella', '#37a9c0'),
    ('Cronobacter', '#35ad9b'),
    ('Psychrobacter', '#3aa5de'),
    ('Cereibacter', '#ce8f31'),
    ('Thiomicrorhabdus', '#a38cf4'),
    ('Pseudogemmobacter', '#f77732'),
    ('Haematospirillum', '#9da131')
]

# Create SVG
svg_lines = ['<svg width="300" height="1200" xmlns="http://www.w3.org/2000/svg">']

for i, (genus, color) in enumerate(genera):
    y = 20 + i * 20
    svg_lines.append(f'  <circle cx="10" cy="{y}" r="8" fill="{color}"/>')
    svg_lines.append(f'  <text x="25" y="{y+4}" font-family="Arial" font-size="12">{genus}</text>')

svg_lines.append('</svg>')

# Write to file
with open('genus_colors.svg', 'w') as f:
    f.write('\n'.join(svg_lines))

print("SVG saved as genus_colors.svg")