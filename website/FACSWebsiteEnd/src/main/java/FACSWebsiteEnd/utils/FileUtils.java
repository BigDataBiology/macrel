package FACSWebsiteEnd.utils;

import FACSWebsiteEnd.Entity.FacsOutIdsTsv;
import FACSWebsiteEnd.Entity.FacsOutTsv;
import FACSWebsiteEnd.Entity.FileInfo;
import FACSWebsiteEnd.common.Constant;
import org.springframework.web.multipart.MultipartFile;

import java.io.IOException;
import java.io.InputStream;
import java.util.*;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

/**
 * @Author: HiramHe
 * @Date: 2019/12/2 17:11
 * QQ:776748935
 */
public class FileUtils {

    public static Map getFileInformation(MultipartFile file){

        Map<String,String> infoMap = new HashMap<String,String>();

        String originalFilename = file.getOriginalFilename();

        // 获取上传文件的类型
        String extension = originalFilename.indexOf(".") != -1 ?
                originalFilename.substring(originalFilename.lastIndexOf(".")+1) : null;

        int lengthOfSuffix = 0;
        if (extension != null){
            lengthOfSuffix = extension.length() + 1;
        }

        String filenameWithOutExtension = originalFilename.substring(0,originalFilename.length()-lengthOfSuffix);

        infoMap.put("filenameWithExtension",originalFilename);
        infoMap.put("filenameWithOutExtension",filenameWithOutExtension);
        infoMap.put("extension",extension);

        return infoMap;
    }

    public static FileInfo Info4TextToFile(String dir,String extension){
        String dot = ".";
        String filenameWithOutExtension = Constant.TEXTFILEPREX + CommonUtils.getUUID();
        String filenameWithExtension =  filenameWithOutExtension + dot + extension;
        String path = dir + filenameWithExtension;

        FileInfo fileInfo = new FileInfo(filenameWithExtension,filenameWithOutExtension, dir, path, extension);
        return fileInfo;
    }

    public static FileInfo setInfo4File(MultipartFile file, String dir){
        FileInfo fileInfo = new FileInfo();

        Map information = getFileInformation(file);

        String filenameWithExtension = information.get("filenameWithExtension").toString();
        String filenameWithOutExtension = information.get("filenameWithOutExtension").toString();
        String extension = information.get("extension").toString();

        // 对文件重命名
        filenameWithOutExtension = filenameWithOutExtension + "-" + CommonUtils.getUUID();
        filenameWithExtension = filenameWithOutExtension + "." + extension;

        String path = dir + filenameWithExtension;

        fileInfo.setFilenameWithExtension(filenameWithExtension);
        fileInfo.setFilenameWithOutExtension(filenameWithOutExtension);
        fileInfo.setDir(dir);
        fileInfo.setPath(path);
        fileInfo.setExtension(extension);

        return fileInfo;
    }

    public static List saveGZInputstreamToObject(InputStream inputStream,Object object){

        List<Object> objects = new ArrayList<>();

        GZIPInputStream gzipInputStream = null;
        try {
            // 解压工作流
            gzipInputStream = new GZIPInputStream(inputStream);
            Scanner scanner = new Scanner(gzipInputStream);

            String line = null;
            String[] words;

            String lineRegex = null;
            Pattern pattern = null;
            Boolean hasHeader = true;
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
            return objects;
        } catch (IOException e) {
            e.printStackTrace();
            return null;
        } finally {
            if (gzipInputStream != null){
                try {
                    gzipInputStream.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

    }


}
