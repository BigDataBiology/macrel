package FACSWebsiteEnd.utils;

import org.springframework.web.multipart.MultipartFile;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @Author: HiramHe
 * @Date: 2019/11/29 11:27
 * QQ:776748935
 */
public class EffectiveCheckUtils {

    public static Boolean strEffectiveCheck(String str){
        return str != null && str.length() >0;
    }

    public static Boolean fileEffectiveCheck(MultipartFile file){
        return file != null && !file.isEmpty();
    }

    public static Boolean emailEffectiveCheck(String email){
        String regEx = "^([a-z0-9A-Z]+[-|\\.]?)+[a-z0-9A-Z]@([a-z0-9A-Z]+(-[a-z0-9A-Z]+)?\\.)+[a-zA-Z]{2,}$";
        Pattern pattern = Pattern.compile(regEx);
        Matcher matcher = pattern.matcher(email);
        if (matcher.matches()){
            return true;
        } else {
            return false;
        }
    }
}
