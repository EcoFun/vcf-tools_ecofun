include("readfasta.jl")

# one line in the VCF file
type SNP
	scaffold
	pos
	id
	ref
	alt
	qual
	kind
	info
	gts
end

# constructor
function SNP(line)
	fields = split(line)

	@assert length(fields)>9 "more than 9 fields in VCF line"

	pos = int(fields[2])
	qual = int(fields[6])
	kind = fields[7]
	gts = fields[10:end]

	SNP(fields[1], pos, fields[3], fields[4], fields[5], qual, kind, fields[8], gts)
end

# scaffold - pos for snp as string
coord(snp) = "$(snp.scaffold) - $(snp.pos)"

isrefcall(snp) = snp.kind == "REFCALL"
callkind(snp) = snp.kind == "REFCALL" ? "REFCALL" : "variant"

validkind(snp) = snp.kind == "REFCALL" || snp.kind == "PASS" || snp.kind == "HapScore"

endpos(snp) = int(split(match(r"END=[0-9]+", snp.info).match, '=')[2])



# merged data
type Merged
	# one block of data from the ref fasta
	chunk
	# keep this one for output
	firstpos
	# last parsed position
	lastpos
	name
	# array of haplotypes
	data
end

# constructor
Merged(chunk) = Merged(chunk, chunk.pos_min, chunk.pos_min-1, chunk.id, Vector{Char}[])


type Config
	minGQ
	minNR
	NACN
	minQUAL
	minRefQUAL
end

# fill up the end with Ns
# !! only call this after you are done with this merged !!
function pad_end(merged)
	missing = merged.chunk.pos_max - merged.lastpos
	
	if missing < 1 
		return
	end

	println("padding end of chunk $(merged.lastpos)-$(merged.chunk.pos_max)")

	to_add = fill('N', missing)

	for ind in merged.data
		append!(ind, to_add)
	end
end

# write data to vector at an offset (resizing if necessary)
function set_at!{T}(vector::Vector{T}, data::Vector{T}, offset = 0, num = 1)
	if offset == 0
		append!(vector, data)
		return
	end

	rng_beg = length(vector) + offset + 1
	rng_end = length(vector) + offset + num
	# splice doesn't like indices out of range
	if rng_end > length(vector)
		rng_end = length(vector)
	end

	@assert rng_end >= rng_beg

	# remove els in range, insert data
	splice!(vector, rng_beg:rng_end, data)
end

# merge one VCF line into the data
function merge(merged::Merged, snp::SNP, config::Config, inds)
	# first line, need to create arrays
	if isempty(merged.data)
		# create arrays for haplotype data
		for i in 1:2*length(snp.gts)
			push!(merged.data, Char[])
		end
	end

	@assert length(merged.data) == 2*length(snp.gts)

# *** global filter

	refcall = isrefcall(snp)
	if snp.qual < (refcall ? config.minRefQUAL : config.minQUAL)
		println(coord(snp) * ": filtered " * 
				"$(callkind(snp))" * " - call quality < " *
				"$(refcall ? config.minRefQUAL : config.minQUAL)")
		# will be filled with 'N' the next time around
		return
	end


	# the actual start of the snp (in ref coordinates)
	snp_begin = snp.pos
	# needed more often, let's just calculate it once
	snp_end = refcall ? endpos(snp) : snp.pos + length(snp.ref) - 1


# *** snp.pos != merged.lastpos+1 (i.e. gaps, overlaps)

	# where the current data starts relative to the end of the *merged* sequence
	offset = 0

	# gap
	if snp_begin-1 > merged.lastpos
		println(coord(snp) * ": filling gap $(merged.lastpos)-$(snp_begin-1)")
		Ns = fill('N', snp_begin-1-merged.lastpos)
		# fill the gap 
		for ind in merged.data
			append!(ind, Ns)
		end
		merged.lastpos = snp_begin-1
	# overlap
	elseif snp_begin ≤ merged.lastpos
		if refcall
			println(coord(snp) * 
				": WARNING, refcall pos smaller than the highest pos processed so far (snp or end pos of a refcall)" *
				" $(merged.lastpos)")
			# refcall fully in range of previous refcall, ignore
			if snp_end <= merged.lastpos
				println(coord(snp) * ": WARNING, redundant refcall (fully in range of previous refcall)")
				return
			end
			# we just start the refcall after the end of the last pos
			# we don't use the offset here so as not to overwrite SNPs
			snp_begin = merged.lastpos + 1
		else
			println(coord(snp) * 
				": WARNING, non-refcall pos smaller than the highest pos processed so far (snp or end pos of a refcall)" *
				" $(merged.lastpos)")

			if length(snp.ref) < length(snp.alt)
				println(coord(snp) * 
					": WARNING, insertion overlapping with refcall")
			end

			# snp starts before current end of sequence
			offset = snp_begin - (merged.lastpos+1)
		end
	end

	@assert merged.lastpos+1 + offset == snp_begin
	@assert snp_end >= snp_begin
	@assert offset <= 0

	# how
	nrepl = snp_end - snp_begin + 1

# *** create ref, alt, adjust merged.lastpos

	ref = Char[]
	alt = Char[]

	# REFCALL
	if refcall
		# get sites from reference fasta
		ref = get_at_pos(merged.chunk, (merged.lastpos+1):snp_end)
		println(coord(snp) * ": REFCALL, $(length(ref)) sites")
	# non-REFCALL
	else
		# convert string into array of chars (for efficiency)
		ref = collect(snp.ref)
		alt = collect(snp.alt)

		diff_l = length(ref) - length(alt)

		# pad ref or alt with '-' if necessary
		if diff_l < 0
			append!(ref, fill('-', -diff_l))
		elseif diff_l > 0
			append!(alt, fill('-', diff_l))
		end
	end

	# only adjust if this is not an overlap
	if snp_end > merged.lastpos 
		merged.lastpos = snp_end
	end

	# create them once for efficiency
	Ns = fill('N', length(ref))
	dashes = fill('-', length(ref))

# *** extend haplotypes according to genotype

	gts = snp.gts

	for ind in 1:length(gts)
		hap = 2 * ind - 1
		
		gt = split(gts[ind], ':')

		# unknown format
		# any number other than 0 and 1 is unexpected for bi-allelic polymorphism
		if !ismatch(r"[0-1][|/][0-1]", gt[1]) && gt[1] != "./." 
			error(coord(snp) * 
				": unexpected genotype $(gt[1]) for individual $ind ($(inds[ind]))")
		end

		GQ = refcall ? 0 : int(gt[4])
		NR = int(gt[5])
		CN = gt[7]		# string is fine for some cases

		# missing data
		if !refcall && gt[1] == "./."
			if CN != "0"
				println(coord(snp) * 
					": WARNING: ./. outside of REFCALL with CN=$(gt[7])" *
					" (GQ=$GQ, NR=$NR) in individual $ind ($(inds[ind])) - put 'N'")
					
					# append Ns
					set_at!(merged.data[hap], Ns, offset, nrepl)
					set_at!(merged.data[hap+1], Ns, offset, nrepl)
			else
				println(coord(snp) * 
					": found ./. outside of REFCALL with CN=0 (GQ=$GQ, NR=$NR)" *
					" in individual $ind ($(inds[ind])) - put '-'")
					
					# append dashes
					set_at!(merged.data[hap], dashes, offset, nrepl)
					set_at!(merged.data[hap+1], dashes, offset, nrepl)
			end
			
			# done with this one
			continue
		end

		# actual, real data

		# do that here, in the other cases we won't need the numerical value
		CN = gt[7] == "NA" ? config.NACN : float(gt[7])

		# filter non-refcalls
		# TODO: check whether any of these need to be tested for REFCALLs
		if !refcall && (GQ < config.minGQ || NR < config.minNR)
			# log message
			if GQ < config.minGQ && NR < config.minNR
				println(coord(snp) * 
					": filtered - GQ<$(config.minGQ) and NR<$(config.minNR) in" *
					" invidual $ind ($(inds[ind]))")
			elseif GQ < config.minGQ
				println(coord(snp) * 
					": filtered - GQ<$(config.minGQ) in" *
					" invidual $ind ($(inds[ind]))")
			elseif NR < config.minNR
				println(coord(snp) * ": filtered - NR<$(config.minNR) in" *
				" invidual $ind ($(inds[ind]))")
			end

			set_at!(merged.data[hap], Ns, offset, nrepl)
			set_at!(merged.data[hap+1], Ns, offset, nrepl)

			continue
		end

		if refcall
			set_at!(merged.data[hap], ref, offset, nrepl)
			set_at!(merged.data[hap+1], ref, offset, nrepl)
		else
			use_ref = int(split(gt[1], ['/', '|']))

			set_at!(merged.data[hap], use_ref[1]==0 ? ref : alt, offset, nrepl)
			set_at!(merged.data[hap+1], use_ref[2]==0 ? ref : alt, offset, nrepl)
		end
	end
end

