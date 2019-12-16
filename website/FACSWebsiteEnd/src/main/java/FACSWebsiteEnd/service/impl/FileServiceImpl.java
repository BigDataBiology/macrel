package FACSWebsiteEnd.service.impl;

import FACSWebsiteEnd.Entity.FacsOutIdsTsv;
import FACSWebsiteEnd.Entity.FileInfo;
import FACSWebsiteEnd.Entity.FacsOutTsv;
import FACSWebsiteEnd.common.Constant;
import FACSWebsiteEnd.service.FileService;
import FACSWebsiteEnd.utils.CommonUtils;
import FACSWebsiteEnd.utils.EffectiveCheckUtils;
import FACSWebsiteEnd.utils.FileUtils;
import org.springframework.stereotype.Service;
import org.springframework.web.multipart.MultipartFile;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

/**
 * @Author: HiramHe
 * @Date: 2019/11/28 16:55
 * QQ:776748935
 */

@Service
public class FileServiceImpl implements FileService {

    @Override
    public FileInfo uploadFileToLocal(MultipartFile file, String savedDir) {

        FileInfo fileInfo = new FileInfo();

        Map information = FileUtils.getFileInformation(file);

        String filenameWithExtension = information.get("filenameWithExtension").toString();
        String filenameWithOutExtension = information.get("filenameWithOutExtension").toString();
        String extension = information.get("extension").toString();

        // 对文件重命名
        filenameWithOutExtension = filenameWithOutExtension + "-" + CommonUtils.getUUID();
        filenameWithExtension = filenameWithOutExtension + "." + extension;

        String path = savedDir + filenameWithExtension;

        fileInfo.setFilenameWithExtension(filenameWithExtension);
        fileInfo.setFilenameWithOutExtension(filenameWithOutExtension);
        fileInfo.setDir(savedDir);
        fileInfo.setPath(path);
        fileInfo.setExtension(extension);

        BufferedOutputStream outputStream = null;
        try {

            outputStream = new BufferedOutputStream(new FileOutputStream(path));
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
    public FileInfo saveTextToFile(String text, String savedDir, String extension) {

        String dot = ".";
        String filenameWithOutExtension = Constant.TEXTFILEPREX + CommonUtils.getUUID();
        String filenameWithExtension =  filenameWithOutExtension + dot + extension;
        String path = savedDir + filenameWithExtension;

        FileInfo fileInfo = new FileInfo(filenameWithExtension,filenameWithOutExtension, savedDir, path, extension);

        File outputFile = new File(path);

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
    public List<Object> readLocalTsvGzToObject(String filePath, Object object){

        // 将读取的tsv的每一行保存到对象中，再将对象放到集合中返回
        List<Object> objects = new ArrayList<Object>();

        // try start.
        try {
            // 文件输入流
            FileInputStream fileInputStream = new FileInputStream(filePath);
            // 解压工作流
            GZIPInputStream gzipInputStream = new GZIPInputStream(fileInputStream);
            Scanner scanner = new Scanner(gzipInputStream);

            String line = null;
            String[] words;

            String lineRegex = null;
            Pattern pattern = null;
            Boolean hasHeader = null;
            if (object instanceof FacsOutTsv){
                lineRegex = "\\s+";
                pattern = Pattern.compile(lineRegex);
                hasHeader = true;
            } else if (object instanceof FacsOutIdsTsv){
                lineRegex = ">|\\s+";
                pattern = Pattern.compile(lineRegex);
                hasHeader = false;
            }

            // while start.
            int i = 0;
            while (scanner.hasNextLine()){
                // 读一行
                line = scanner.nextLine();
                i++;

                if (i == 1){
                    if (hasHeader){
                        continue;
                    }
                }

//                if (i>3){
//                    break;
//                }

                if (object instanceof FacsOutTsv){
                    // 利用正则表达式，将该行分割成各个单词
                    words = pattern.split(line);
                    if (EffectiveCheckUtils.arrayEffectiveCheck(words)){
                        String access = words[0];
                        String sequence = words[1];
                        String amp_family = words[2];
                        Double amp_probability = Double.valueOf(words[3]);
                        String hemolytic = words[4];
                        Double hemolytic_probability = Double.valueOf(words[5]);

                        FacsOutTsv facsOutTsv = new FacsOutTsv(access,sequence,amp_family,amp_probability,hemolytic,hemolytic_probability);
                        objects.add(facsOutTsv);

                    }
                } else if (object instanceof FacsOutIdsTsv){
                    words = pattern.split(line);
                    if (EffectiveCheckUtils.arrayEffectiveCheck(words)){
                        // todo 正则表达式分割后有一个空值
                        String smORF_num = words[1];
                        String sequence = words[2];
                        String uniqueMark = words[3];

                        FacsOutIdsTsv facsOutIdsTsv = new FacsOutIdsTsv(smORF_num,sequence,uniqueMark);
                        objects.add(facsOutIdsTsv);

                    }
                }
            }
            // while end.

        } catch (Exception e){
            e.printStackTrace();
        }
        // try end.

        return objects;
    }
}
