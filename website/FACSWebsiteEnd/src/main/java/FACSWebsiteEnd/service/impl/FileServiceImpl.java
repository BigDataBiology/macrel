package FACSWebsiteEnd.service.impl;

import FACSWebsiteEnd.Entity.FileInfo;
import FACSWebsiteEnd.config.RemoteProperties;
import FACSWebsiteEnd.service.FileService;
import FACSWebsiteEnd.utils.CommandUtils;
import FACSWebsiteEnd.utils.FileUtils;
import org.springframework.stereotype.Service;
import org.springframework.web.multipart.MultipartFile;

import java.io.*;
import java.util.List;

/**
 * @Author: HiramHe
 * @Date: 2019/11/28 16:55
 * QQ:776748935
 */

@Service
public class FileServiceImpl implements FileService {

    @Override
    public FileInfo uploadFileToLocal(MultipartFile file, String savedDir) {

        FileInfo fileInfo = FileUtils.setInfo4File(file,savedDir);

        BufferedOutputStream outputStream = null;
        try {

            outputStream = new BufferedOutputStream(new FileOutputStream(fileInfo.getPath()));
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
    public FileInfo saveTextToFileLocally(String text, String savedDir, String extension) {

        FileInfo fileInfo = FileUtils.Info4TextToFile(savedDir, extension);

        File outputFile = new File(fileInfo.getPath());

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

    @Override
    public List readFileToObjectFromLocal(String filePath, Object object){

        // 将读取的tsv的每一行保存到对象中，再将对象放到集合中返回
        List objects = null;

        try {
            // 文件输入流
            FileInputStream fileInputStream = new FileInputStream(filePath);
            objects = FileUtils.saveGZInputstreamToObject(fileInputStream,object);

        } catch (Exception e){
            e.printStackTrace();
        }

        return objects;
    }

    @Override
    public List readFileToObjectFromRemote(RemoteProperties remoteProperties, String filePath, Object object) {
        List objects;
        objects = CommandUtils.downloadFileFromRemote(remoteProperties, filePath, object);

        return objects;
    }

    @Override
    public FileInfo uploadFileToRemote(RemoteProperties remoteProperties, String dir, MultipartFile file) {

        FileInfo fileInfo = FileUtils.setInfo4File(file,dir);
        Boolean isSuccessfull = CommandUtils.uploadFileToRemote(remoteProperties, fileInfo.getDir(), fileInfo.getFilenameWithExtension(), file);
        return fileInfo;
    }

    @Override
    public FileInfo saveContentToFileRemotely(RemoteProperties remoteProperties, String dir, String extension, String content) {

        FileInfo fileInfo = FileUtils.Info4TextToFile(dir,extension);
        boolean isSuccessfull = CommandUtils.saveContentToFileRemotely(remoteProperties,fileInfo.getDir(),fileInfo.getFilenameWithExtension(),content);
        if (isSuccessfull){
            return fileInfo;
        } else {
            return null;
        }
    }
}
