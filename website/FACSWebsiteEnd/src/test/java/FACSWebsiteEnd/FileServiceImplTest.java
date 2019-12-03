package FACSWebsiteEnd;

import FACSWebsiteEnd.service.FileService;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.test.context.junit4.SpringRunner;

/**
 * @Author: HiramHe
 * @Date: 2019/11/29 12:36
 * QQ:776748935
 */

@RunWith(SpringRunner.class)
@SpringBootTest
public class FileServiceImplTest {

    @Autowired
    FileService fileService;

    @Test
    public void testSaveTextToFile(){
        String text = "I am hhm," +
                "who are you";
        String extension = "fastq";

        fileService.saveTextToFile(text,extension);
    }
}
