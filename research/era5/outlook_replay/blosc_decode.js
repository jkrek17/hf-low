// Decode a Blosc1 frame (zstd, typesize 1, bitshuffle flag handled in Python) from argv[2] to raw bytes at argv[3].
const fs=require('fs'),zlib=require('zlib');
const b=fs.readFileSync(process.argv[2]);
const flags=b[2],typesize=b[3],nbytes=b.readUInt32LE(4),bs=b.readUInt32LE(8);
const comp=(flags>>5)&7, nblocks=Math.ceil(nbytes/bs);
if(comp!==4||typesize!==1) throw new Error('unsupported frame '+flags+' '+typesize);
const out=Buffer.alloc(nbytes);
for(let i=0;i<nblocks;i++){
  let pos=b.readUInt32LE(16+4*i);
  const bsz=Math.min(bs,nbytes-i*bs);
  const nsplit=(flags&0x10)?1:((bsz>=128)?1:1);
  const cs=b.readInt32LE(pos); pos+=4;
  if(cs===bsz){b.copy(out,i*bs,pos,pos+cs)}
  else{const d=zlib.zstdDecompressSync(b.subarray(pos,pos+cs));if(d.length!==bsz)throw new Error('size');d.copy(out,i*bs)}
}
fs.writeFileSync(process.argv[3],out);
