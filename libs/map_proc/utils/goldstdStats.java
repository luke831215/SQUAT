import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Date;

public class goldstdStats {
	public static int[] parse(final String answer) {
		// length, #errors, max cont. #errors
		int[] rst = new int[] { 0, 0, Integer.MIN_VALUE };

		int cont_cnt = 0;

		StringBuilder pos_in_char = new StringBuilder();

		for (int i = 0; i < answer.length() + 1; i++) {
			char tmp = 'X';

			if (i < answer.length()) {
				tmp = answer.charAt(i);
			}

			if (Character.isDigit(tmp)) {
				pos_in_char.append(tmp);
			} else {
				int dist = Integer.valueOf(pos_in_char.toString());

				rst[0] += dist;
				if (i < answer.length()) {
					// the symbol itself
					rst[0] += 1;
					rst[1]++;
				}

				if (dist == 0) {
					cont_cnt++;
				} else {
					if (cont_cnt > 0) {
						rst[2] = Math.max(rst[2], cont_cnt);
					}

					cont_cnt = 1;
				}

				pos_in_char.delete(0, pos_in_char.length());
			}
		}

		if (rst[2] < 0) {
			rst[2] = 0;
		}

		return rst;
	}

	public static void main(String[] args) throws IOException {
		if (args.length < 1) {
			System.err.println("Using: " + goldstdStats.class.getName()
					+ " [Gold standard (ans)]");
			return;
		}

		long total_bases = 0;
		int total_error = 0;

		int cnt_read_has_error = 0;
		int cnt_read_no_error = 0;

		int error_min = Integer.MAX_VALUE;
		int error_max = Integer.MIN_VALUE;

		int cont_error_max = Integer.MIN_VALUE;

		BufferedWriter br_out = new BufferedWriter(
				new FileWriter(args[0] + ".report"));

		// for report output
		br_out.write("Statistic of " + args[0] + " on "
				+ (new SimpleDateFormat("yyyy-MM-dd HH:mm")).format(new Date())
				+ "\n");
		br_out.write("\n");

		BufferedReader br_ans = new BufferedReader(new FileReader(args[0]));
		while (br_ans.ready()) {
			String ans = br_ans.readLine().trim();
			if (!ans.equals("")) {
				String[] data = ans.split("\t", 2);

				int[] rst = parse(data[1]);

				total_bases += rst[0];
				total_error += rst[1];

				if (rst[1] == 0) {
					cnt_read_no_error++;
				} else {
					cnt_read_has_error++;
				}

				// min and max
				error_min = Math.min(error_min, rst[1]);
				error_max = Math.max(error_max, rst[1]);

				cont_error_max = Math.max(cont_error_max, rst[2]);

				// for report output
				br_out.write(String.format("%10s", data[0]));
				for (int i = 0; i < rst.length; i++) {
					br_out.write(String.format("\t%4s", rst[i]));
				}
				br_out.write("\t" + data[1] + "\n");
			}
		}
		br_ans.close();

		br_out.close();

		// for stdout
		System.out.print(total_bases + "\t" + total_error);
		System.out.print("\t" + cnt_read_no_error + "\t" + cnt_read_has_error);
		System.out.print("\t" + error_min + "\t" + error_max);
		System.out.print("\t" + cont_error_max);
		System.out.println();
	}
}

