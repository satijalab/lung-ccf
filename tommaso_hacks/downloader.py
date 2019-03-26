""""
Downloader from Box folder w/o Box client

- Install boxsdk from [github](https://github.com/box/box-python-sdk)
- Create custom app on [website](https://developer.box.com/) and use recommended authentication system.
- Set up JWT followimg instruction on [webpage](https://developer.box.com/docs/setting-up-a-jwt-app).
- Add enterprise access to box account.
"""

import os
import logging
import boxsdk as box



def download_items(folder_content, folder_name='./HuBMAP_Satija_Kharchenko_CCF'):
    """
        Recursively create folders and download files.
    """
    print('Entering {}.'.format(folder_name))
    for item in folder_content:
        if hasattr(item, 'download_to') is True:
            file_name = os.path.join(folder_name, item.name)
            with open(file_name, 'bw') as file:
                item.download_to(file)
                logging.info('Downloaded {}.'.format(file))
        else:
            new_folder_name = os.path.join(folder_name, item.name)
            os.makedirs(new_folder_name)
            try:
                new_folder_content = item.get_items()
                download_items(new_folder_content, folder_name=new_folder_name)
            except:
                print('Error on {}'.format(item))


if __name__ == '__main__':

    # Authenticate
    sdk = box.JWTAuth.from_settings_file('/home/tbiancal/.ssh/hubmap-downloader_box.json')
    client = box.Client(sdk)

    # Get root folder
    root_folder = client.folder(folder_id='0')
    root_items = root_folder.get_items()
    items = [item for item in root_items]
    folder = items[0]
    assert folder.name == 'HuBMAP_Satija_Kharchenko_CCF'
    items = folder.get_items()

    download_items(items)
