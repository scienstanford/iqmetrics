function [dinfo, failure, err_rate] = ecd_info(cinfo)

% BCH error correction
[dinfo err] = decode(cinfo, 15, 5, 'bch');
err_rate = sum(err(1:5:end))/length(err);

% Cyclic Redundancy Check (CRC)
dcrc = dinfo(163:178);
dinfo = dinfo(1:162)';
dinfoa = [dinfo 0 0 0 0 0 0];
fid = fopen('CRC\dinfo.dat', 'wb');
fwrite(fid, uint8(dinfoa), 'uint8');
fclose(fid);
command = ['CRC\CRC CRC\dinfo.dat CRC\ncrc.dat'];
dos(command);
fid = fopen('CRC\ncrc.dat', 'rb');
crc = fread(fid);
fclose(fid);

failure = sum(abs(crc-dcrc)) > 0;

return