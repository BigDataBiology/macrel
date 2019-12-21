package FACSWebsiteEnd.service.impl;

import FACSWebsiteEnd.Entity.FileInfo;
import FACSWebsiteEnd.common.ResultCode;
import FACSWebsiteEnd.common.ResultObject;
import FACSWebsiteEnd.config.RemoteProperties;
import FACSWebsiteEnd.service.FileService;
import FACSWebsiteEnd.utils.CommandUtils;
import FACSWebsiteEnd.utils.FileUtils;
import org.omg.PortableInterceptor.INACTIVE;
import org.springframework.stereotype.Service;
import org.springframework.web.multipart.MultipartFile;

import javax.servlet.http.HttpServletResponse;
import java.io.*;
import java.net.URLEncoder;
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
        objects = CommandUtils.downloadFileToObjectFromRemote(remoteProperties, filePath, object);

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

    @Override
    public ResultObject getFileForDownloadFromLocal(File file, String filename, HttpServletResponse response) {

        FileInputStream fileInputStream = null;
        try {
            fileInputStream = new FileInputStream(file);
            Boolean isSuccessful = writeFileToBrowser(response,fileInputStream,filename);

            if (isSuccessful){
                return null;
            } else {
                System.out.println("fail to write file to browser.");
                return ResultObject.failure(ResultCode.SERVER_INTERNAL_ERROR);
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
            return ResultObject.failure(ResultCode.FILE_NOT_EXIST);
        } finally {
            if (fileInputStream != null){
                try {
                    fileInputStream.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

    }

    @Override
    public ResultObject getFileForDownloadFromRemote(RemoteProperties remoteProperties, String filePath, String filename, HttpServletResponse response) {

        InputStream inputStream = null;
        try {

            inputStream = CommandUtils.getFileForDownloadFromRemote(remoteProperties, filePath);

            if (inputStream != null){
                Boolean isSuccessful = writeFileToBrowser(response,inputStream,filename);
                if (isSuccessful){
                    return null;
                } else {
                    System.out.println("fail to write file to browser.");
                    return ResultObject.failure(ResultCode.SERVER_INTERNAL_ERROR);
                }
            } else {
                System.out.println("未能读取远程服务器上的文件.");
                return ResultObject.failure(ResultCode.FILE_NOT_EXIST);
            }

        } finally {
            if (inputStream != null){
                try {
                    inputStream.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            CommandUtils.disConnectChannel();
        }

    }

    private Boolean writeFileToBrowser(HttpServletResponse response, InputStream inputStream, String filename){

        // 实现文件下载
        byte[] buffer = new byte[1024];
        BufferedInputStream bufferedInputStream = null;
        OutputStream outputStream = null;
        try {

            //1.设置文件ContentType类型，这样设置，会自动判断下载文件类型
            response.setContentType("multipart/form-data");
            //2.设置文件头：最后一个参数是设置下载文件名
            String newFilename = URLEncoder.encode(filename, "UTF-8");
            response.setHeader("Content-Disposition", "attachment;filename=" + newFilename);
            // 便于前端获取文件名
            response.setHeader("fileName", newFilename);
            response.setHeader("Access-Control-Expose-Headers", "fileName");

            bufferedInputStream = new BufferedInputStream(inputStream);

            //3.通过response获取OutputStream对象
            outputStream = response.getOutputStream();

            int i = bufferedInputStream.read(buffer);
            while (i != -1) {
                outputStream.write(buffer, 0, i);
                i = bufferedInputStream.read(buffer);
            }

            return true;

        } catch (IOException e) {
            e.printStackTrace();
            return false;
        } finally {
            if (bufferedInputStream != null){
                try {
                    bufferedInputStream.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }
}
