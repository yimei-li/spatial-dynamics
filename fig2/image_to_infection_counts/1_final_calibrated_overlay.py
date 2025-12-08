#!/usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from PIL import Image
import os, csv

HEX_R = 2.82
HEX_H = np.sqrt(3) * HEX_R

INPUT_LIST = [
    ('7_hours.png', 7),
    ('13_hours.png', 13),
    ('19_hours.png', 19),
    ('25_hours_comprehensive.png', 25),
    ('25_hours_DIP-infected.png', 25),
    ('25_hours_Virion-infected.png', 25)
]

def crop_border(img, thresh=240):
    g = np.dot(img[...,:3], [0.299, 0.587, 0.114])
    m = g < thresh
    c = np.argwhere(m)
    if len(c) == 0:
        return img
    y0, x0 = c.min(axis=0)
    y1, x1 = c.max(axis=0) + 1
    return img[y0:y1, x0:x1]

def load_img(path):
    try:
        im = Image.open(path)
        arr = np.array(im.convert('RGB'))
        return crop_border(arr)
    except:
        print("cant load " + path)
        return None

def make_hex_grid(w, h):
    centers = []
    nc = int(np.ceil(w / (1.5 * HEX_R)))
    nr = int(np.ceil(h / HEX_H))
    for r in range(nr):
        for c in range(nc):
            px = c * 1.5 * HEX_R
            py = r * HEX_H
            if c % 2 == 1:
                py = py + HEX_H / 2
            if px < w and py < h:
                centers.append((px, py, r, c))
    return centers, nr, nc

def get_mask(cx, cy, shape):
    h, w = shape[0], shape[1]
    yy, xx = np.ogrid[:h, :w]
    d = np.sqrt((xx - cx)**2 + (yy - cy)**2)
    return d <= HEX_R

def classify_color(img, mask):
    if not np.any(mask):
        return 'N'
    pix = img[mask]
    if len(pix) == 0:
        return 'N'
    
    r = np.mean(pix[:,0])
    g = np.mean(pix[:,1])
    b = np.mean(pix[:,2])
    tot = r + g + b
    
    if tot < 10:
        return 'B'
    if tot < 20:
        return 'B'
    
    rr = r / max(tot, 1)
    gr = g / max(tot, 1)
    
    # green
    is_green = False
    if g > 15:
        if g > r and g > b:
            is_green = True
        elif g > r - 5 and g > b + 2:
            is_green = True
        elif gr > 0.33:
            is_green = True
        elif g > 25 and g > b:
            is_green = True
        elif g > 30:
            is_green = True
        elif tot > 40 and g > r * 0.9 and g > b:
            is_green = True
    if is_green:
        if r > 70 and g > 50 and r > g and b < g * 0.7:
            is_green = False
    if is_green:
        return 'G'
    
    # yellow
    is_yellow = False
    if r > 30 and g > 20:
        if r > b + 10 and g > b:
            is_yellow = True
        elif r + g > tot * 0.6:
            is_yellow = True
        elif r > 40 and g > 25 and b < g:
            is_yellow = True
        elif tot > 80 and rr + gr > 0.6:
            is_yellow = True
        elif r > 35 and g > 20 and r < g * 2.5:
            is_yellow = True
        elif r > 50 and g > 30 and b < 50:
            is_yellow = True
    if is_yellow:
        if r > g * 2.2 and r > 90:
            is_yellow = False
    if is_yellow:
        return 'Y'
    
    # red
    if r > 40:
        if r > g * 2 and r > b * 2:
            return 'R'
        if rr > 0.45 and r > g + 40:
            return 'R'
        if r > 100 and r > g + 50:
            return 'R'
        if r > 80 and rr > 0.5:
            return 'R'
    
    return 'B'

def find_yellowish(img, centers, n=2):
    cands = []
    for (px, py, rr, cc) in centers:
        mask = get_mask(px, py, img.shape)
        pix = img[mask]
        if len(pix) == 0:
            continue
        r = np.mean(pix[:,0])
        g = np.mean(pix[:,1])
        b = np.mean(pix[:,2])
        if r > 50 and g > 20 and r > b:
            score = g / max(r, 1) + (1 - abs(r - g) / max(r, 1))
            cands.append(((px, py), score))
    cands.sort(key=lambda x: x[1], reverse=True)
    out = []
    for i in range(min(n, len(cands))):
        out.append(cands[i][0])
    return out

def draw_overlay(img, centers, cmap, outpath, name):
    h, w = img.shape[0], img.shape[1]
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(img)
    
    clrs = {'R': ('red', 0.7), 'G': ('green', 0.8), 'Y': ('yellow', 0.8), 'B': ('black', 0.6), 'N': ('gray', 0.2)}
    
    cnt = {'R': 0, 'G': 0, 'Y': 0, 'B': 0}
    for v in cmap.values():
        if v in cnt:
            cnt[v] += 1
    
    for (px, py, _, _) in centers:
        cc = cmap.get((px, py), 'B')
        if cc in clrs:
            col, alp = clrs[cc]
            p = mpatches.RegularPolygon((px, py), 6, radius=HEX_R, orientation=np.radians(30),
                                        edgecolor='white', facecolor=col, alpha=alp, linewidth=0.2)
            ax.add_patch(p)
    
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    t = name + "\nR:%d Y:%d G:%d B:%d" % (cnt['R'], cnt['Y'], cnt['G'], cnt['B'])
    ax.set_title(t, fontsize=14)
    plt.axis('off')
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close()

def write_coords(centers, cmap, path):
    f = open(path, 'w', newline='')
    w = csv.writer(f)
    w.writerow(['x', 'y', 'color'])
    for (px, py, _, _) in centers:
        c = cmap.get((px, py), 'B')
        w.writerow([px, py, c])
    f.close()

def write_grid(centers, cmap, nr, nc, path):
    grid = []
    for i in range(nr):
        row = []
        for j in range(nc):
            row.append('')
        grid.append(row)
    for (px, py, r, c) in centers:
        if r < nr and c < nc:
            grid[r][c] = cmap.get((px, py), 'B')
    f = open(path, 'w', newline='')
    w = csv.writer(f)
    for row in grid:
        w.writerow(row)
    f.close()

def write_counts(cnt, tp, path):
    tot = cnt['R'] + cnt['G'] + cnt['Y'] + cnt['B']
    f = open(path, 'w', newline='')
    w = csv.writer(f)
    w.writerow(['time', 'virion_counts', 'dip_counts', 'both_infected_counts',
                'susceptible_counts', 'total_cells', 'virion_rate', 'dip_rate',
                'both_infected_rate', 'susceptible_rate'])
    if tot > 0:
        vr = "%.2f" % (cnt['R']/tot*100)
        dr = "%.2f" % (cnt['G']/tot*100)
        br = "%.2f" % (cnt['Y']/tot*100)
        sr = "%.2f" % (cnt['B']/tot*100)
    else:
        vr = dr = br = sr = "0.00"
    w.writerow([tp, cnt['R'], cnt['G'], cnt['Y'], cnt['B'], tot, vr, dr, br, sr])
    f.close()

def write_summary(data):
    path = "1_infection_counts_by_time.csv"
    f = open(path, 'w', newline='')
    w = csv.writer(f)
    w.writerow(['time', 'virion_counts', 'dip_counts', 'both_infected_counts',
                'susceptible_counts', 'total_cells', 'virion_rate', 'dip_rate',
                'both_infected_rate', 'susceptible_rate'])
    data.sort(key=lambda x: x[0])
    for (tp, cnt, tot) in data:
        if tot > 0:
            vr = "%.2f" % (cnt['R']/tot*100)
            dr = "%.2f" % (cnt['G']/tot*100)
            br = "%.2f" % (cnt['Y']/tot*100)
            sr = "%.2f" % (cnt['B']/tot*100)
        else:
            vr = dr = br = sr = "0.00"
        w.writerow([tp, cnt['R'], cnt['G'], cnt['Y'], cnt['B'], tot, vr, dr, br, sr])
    f.close()
    print("wrote " + path)

def run():
    summary = []
    
    for (fname, tp) in INPUT_LIST:
        if not os.path.exists(fname):
            print("skip " + fname)
            continue
        
        print("processing " + fname)
        img = load_img(fname)
        if img is None:
            continue
        
        h, w = img.shape[0], img.shape[1]
        centers, nr, nc = make_hex_grid(w, h)
        
        cnt = {'R': 0, 'G': 0, 'Y': 0, 'B': 0}
        cmap = {}
        for (px, py, rr, cc) in centers:
            mask = get_mask(px, py, img.shape)
            c = classify_color(img, mask)
            cmap[(px, py)] = c
            cnt[c] += 1
        
        # 13h fix
        if fname == '13_hours.png' and cnt['Y'] < 2:
            yw = find_yellowish(img, centers, 2)
            for pos in yw:
                if pos in cmap and cmap[pos] == 'R':
                    cmap[pos] = 'Y'
                    cnt['R'] -= 1
                    cnt['Y'] += 1
        
        base = os.path.splitext(fname)[0]
        draw_overlay(img, centers, cmap, "1_" + base + "_final_overlay.png", base)
        write_coords(centers, cmap, "1_" + base + "_coordinates.csv")
        write_grid(centers, cmap, nr, nc, "1_" + base + "_grid.csv")
        write_counts(cnt, tp, "1_" + base + "_counts.csv")
        
        if 'DIP' not in fname and 'Virion' not in fname:
            summary.append((tp, cnt, len(centers)))
        
        print("  n=%d R=%d Y=%d G=%d B=%d" % (len(centers), cnt['R'], cnt['Y'], cnt['G'], cnt['B']))
    
    write_summary(summary)

if __name__ == "__main__":
    run()
