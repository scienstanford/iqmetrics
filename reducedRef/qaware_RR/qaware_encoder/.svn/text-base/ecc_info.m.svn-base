function einfo = ecc_info(info)

% add Cyclic Redundancy Check (CRC) bits
infoa = [info 0 0 0 0 0 0];
fid = fopen('CRC\info.dat', 'wb');
fwrite(fid, uint8(infoa), 'uint8');
fclose(fid);
command = ['CRC\CRC CRC\info.dat CRC\crc.dat'];
dos(command);
fid = fopen('CRC\crc.dat', 'rb');
crc = fread(fid);
fclose(fid);
cinfo = [info crc'];

% BCH error correction coding
einfo = encode(cinfo, 15, 5, 'bch');

return