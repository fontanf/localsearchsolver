import gdown
import os
import pathlib
import time


def download(id):
    for _ in range(3):
        try:
            gdown.download(id=id, output="data.7z")
            os.system("7z x data.7z -odata")
            pathlib.Path("data.7z").unlink()
        except:
            time.sleep(10)
            continue
        return
    sys.exit(1)


download("1_sbLgEoWCFo4gEsg9BE8eG4De4f8fkF0")
download("1IPmVu2rtrLcDyaLzwml5aiAKLPEylMc8")
download("1PT4Lw48TRSeDEPacHNT8UmHL6o0i_eqN")
download("1grgLEzCJqyWxv1Xb5n1wUJFr7AtlPNQR")
download("1yxAhYw-ViJnJTKVnvjtYOA78ZXQ_T6io")
download("192DpRsSKKhU6UZd8GVdU0ouZKEKfnEcp")
download("1GWDGAQUj3PeO395drzgdY4AcUvy7y75V")
download("1o2cyipclS9sPRJ0Kt7nJBz1yJnDLLOju")
download("1YAeuhKJcF0uOkdqUMASEPv6zONSftAWI")
download("1XHoPvNPg0gJLxkdFrh9xnJqFoPBKZsyv")
download("18kdpQrQQZ1EMRS6MQhz7XZmlH61ECtDs")
download("1ckAQm2PPt3HBvsSPZjiSBta689GMjCzk")
download("12tY8xeFSXs6ybNo-QethEpE_Ue2dkJ6g")
download("14GtIRtWlwzsg3gS7ODLR1EZHD0TH-mtR")
download("1CdYNxVGoh8U4936U5lzpvJqWVRqr7wzg")
