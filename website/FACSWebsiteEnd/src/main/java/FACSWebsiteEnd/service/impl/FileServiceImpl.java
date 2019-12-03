package FACSWebsiteEnd.service.impl;

import FACSWebsiteEnd.Entity.FileInfo;
import FACSWebsiteEnd.common.Constant;
import FACSWebsiteEnd.service.FileService;
import FACSWebsiteEnd.utils.CommonUtils;
import FACSWebsiteEnd.utils.FileUtils;
import org.springframework.stereotype.Service;
import org.springframework.web.multipart.MultipartFile;

import java.io.*;
import java.util.Map;

/**
 * @Author: HiramHe
 * @Date: 2019/11/28 16:55
 * QQ:776748935
 */

@Service
public class FileServiceImpl implements FileService {

    @Override
    public FileInfo upload(MultipartFile file) {

        FileInfo fileInfo = new FileInfo();

        Map information = FileUtils.getFileInformation(file);
        String filename = information.get("filename").toString();
        String type = information.get("type").toString();

        fileInfo.setFilename(filename);
        fileInfo.setExtension(type);

        String fullpath = Constant.FILESAVED_DIR+filename;
        fileInfo.setPath(Constant.FILESAVED_DIR);
        fileInfo.setFullpath(fullpath);

        BufferedOutputStream outputStream = null;
        try {

            outputStream = new BufferedOutputStream(new FileOutputStream(fullpath));
            outputStream.write(file.getBytes());
            outputStream.flush();

            return fileInfo;
        } catch (IOException e) {
            e.printStackTrace();
            return  null;
        } finally {
            try {
                if (outputStream != null) {
                    outputStream.close();
                }
            } catch (IOException e) {
                e.printStackTrace();
            }

        }
    }

    @Override
    public FileInfo saveTextToFile(String text, String extension) {

        String dot = ".";
        String fileName = Constant.FASTQPREX + CommonUtils.getUUID() + dot + extension;
        String fullpath = Constant.FILESAVED_DIR + fileName;

        FileInfo fileInfo = new FileInfo(fileName, Constant.FILESAVED_DIR, fullpath, extension);

        File outputFile = new File(fullpath);

        BufferedWriter bufferedWriter = null;

        try {
            FileWriter fileWriter = new FileWriter(outputFile);
            bufferedWriter = new BufferedWriter(fileWriter);

            bufferedWriter.write(text);
            bufferedWriter.flush();

            return fileInfo;

        } catch (IOException e) {
            e.printStackTrace();
            return null;
        } finally {
            try {
                if (bufferedWriter != null) {
                    bufferedWriter.close();
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
}
