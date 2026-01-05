import os
import re
from reportlab.lib.pagesizes import A4, landscape
from reportlab.pdfgen import canvas
from reportlab.lib.units import cm
from PIL import Image

IMAGE_DIR = "histograms"
OUTPUT_PDF = "histogram_report_landscape_square.pdf"

pattern = re.compile(r"index_(\d+)_silence_(\d+)")

def sort_key(filename):
    m = pattern.search(filename)
    if m:
        return (int(m.group(1)), int(m.group(2)))
    return (999999, 999999)

images = sorted(
    [f for f in os.listdir(IMAGE_DIR) if f.lower().endswith(".png")],
    key=sort_key
)

page_width, page_height = landscape(A4)

# 🔧 DENGELİ BOŞLUKLAR
margin_x = 1.2 * cm
margin_y = 1.2 * cm
inner_padding = 0.3 * cm

cols = 4
rows = 2

cell_width = (page_width - 2 * margin_x) / cols
cell_height = (page_height - 2 * margin_y) / rows

# 🔑 Kareye zorla yaklaştır
cell_size = min(cell_width, cell_height)

c = canvas.Canvas(OUTPUT_PDF, pagesize=landscape(A4))

for i, img_name in enumerate(images):
    pos = i % 8
    col = pos % cols
    row = pos // cols

    x = margin_x + col * cell_width + (cell_width - cell_size) / 2
    y = page_height - margin_y - (row + 1) * cell_height + (cell_height - cell_size) / 2

    img_path = os.path.join(IMAGE_DIR, img_name)
    img = Image.open(img_path)
    iw, ih = img.size
    aspect = iw / ih

    avail = cell_size - 2 * inner_padding

    # Aspect korunur (fit)
    if aspect > 1:
        draw_w = avail
        draw_h = avail / aspect
    else:
        draw_h = avail
        draw_w = avail * aspect

    offset_x = x + (cell_size - draw_w) / 2
    offset_y = y + (cell_size - draw_h) / 2

    c.drawImage(
        img_path,
        offset_x,
        offset_y,
        draw_w,
        draw_h,
        preserveAspectRatio=True
    )

    if pos == 7:
        c.showPage()

if len(images) % 8 != 0:
    c.showPage()

c.save()

print("A4 yatay, kare ve okunaklı PDF oluşturuldu:", OUTPUT_PDF)
