package FACSWebsiteEnd;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.test.context.junit4.SpringRunner;

/**
 * @Author: HiramHe
 * @Date: 2019/12/10 20:52
 * QQ:776748935
 */

public class FileTest {

    @Test
    public void getFilenameFromUriTest(){

        String uri = "/home/HiramHe/tmp/facs_out/excontigs-f64253a61a7b4ba2a0193c41cf524cf7/FACS_OUT.tsv.gz";
        String uri2 = "FACS_OUT.tsv.gz";

        String filenameWithExtension = uri2.lastIndexOf("/") != -1 ?
                uri2.substring(uri2.lastIndexOf("/")+1) : uri2;

        System.out.println(filenameWithExtension);
    }


}
