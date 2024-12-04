import os
import shutil

def draftToMain(window, widget, basePath, dataPath, draftPath):
    for root, dirs, files in os.walk(draftPath):
        for name in files:
            src = os.path.join(root, name)
            dst = os.path.join(root.replace(draftPath, dataPath), name)
            shutil.move(src, dst)
            widget.insert('end', '{src} renombrado a {dst}.\n'.format(src=src, dst=dst))
            widget.see('end')
            window.update_idletasks()

    widget.insert('end', '\nProceso terminado.\n')
    widget.see('end')
    window.update_idletasks()
